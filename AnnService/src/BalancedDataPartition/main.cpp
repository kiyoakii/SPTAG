#include <mpi.h>
#include "inc/Core/Common/DistanceUtils.h"
#include "inc/Core/Common/Dataset.h"
#include "inc/Core/Common/BKTree.h"
#include "inc/Helper/VectorSetReader.h"
#include "inc/Helper/CommonHelper.h"

using namespace SPTAG;

class PartitionOptions : public Helper::ReaderOptions
{
public:
    PartitionOptions():Helper::ReaderOptions(VectorValueType::Float, 0, "|", 32)
    {
        AddRequiredOption(m_inputFiles, "-i", "--input", "Input raw data.");
        AddRequiredOption(m_clusterNum, "-c", "--numclusters", "Number of clusters.");
        AddOptionalOption(m_stopDifference, "-d", "--diff", "Clustering stop center difference.");
        AddOptionalOption(m_maxIter, "-r", "--iters", "Max clustering iterations.");
        AddOptionalOption(m_localSamples, "-s", "--samples", "Number of samples for fast clustering.");
        AddOptionalOption(m_lambda, "-l", "--lambda", "lambda for balanced level.");
        AddOptionalOption(m_distMethod, "-m", "--dist", "Distance method (L2 or Cosine).");
		AddOptionalOption(m_outdir, "-o", "--outdir", "Output directory.");
		AddOptionalOption(m_seed, "-e", "--seed", "Random seed.");
		AddOptionalOption(m_initIter, "-x", "--init", "Number of iterations for initialization.");
    }

    ~PartitionOptions() {}

    std::string m_inputFiles;
    int m_clusterNum;
    
    float m_stopDifference = 0.000001f;
    int m_maxIter = 100;
    int m_localSamples = 1000;
    float m_lambda = 0.000001f;
	int m_seed = -1;
	int m_initIter = 3;
    DistCalcMethod m_distMethod = DistCalcMethod::L2;

    std::string m_centers = "centers.bin";
	std::string m_outdir = "-";
    std::string m_outfile = "vectors.bin";
    std::string m_outmetafile = "meta.bin";
    std::string m_outmetaindexfile = "metaindex.bin";
} options;

template <typename T>
bool LoadCenters(T* centers, SizeType row, DimensionType col) {
    if (fileexists(options.m_centers.c_str())) {
        std::ifstream inputStream(options.m_centers, std::ifstream::binary);
        if (!inputStream.is_open()) {
            fprintf(stderr, "Failed to read center file %s.\n", options.m_centers.c_str());
            return false;
        }

        SizeType r;
        DimensionType c;
        inputStream.read((char*)&r, sizeof(SizeType));
        inputStream.read((char*)&c, sizeof(DimensionType));
        if (r != row || c != col) return false;
        
        inputStream.read((char*)centers, sizeof(T)*row*col);
        inputStream.close();
        std::cout << "load centers(" << row << "," << col << ") from file " << options.m_centers << std::endl << std::flush;
        return true;
    }
    return false;
}

template <typename T>
bool SaveCenters(T* centers, SizeType row, DimensionType col) {
    std::ofstream outputStream(options.m_centers, std::ofstream::binary);
    if (!outputStream.is_open()) {
        fprintf(stderr, "Failed to open center file %s to write.\n", options.m_centers.c_str());
        return false;
    }
    outputStream.write((char*)&row, sizeof(SizeType));
    outputStream.write((char*)&col, sizeof(DimensionType));
    outputStream.write((char*)centers, sizeof(T)*row*col);
    outputStream.close();
    std::cout << "save centers(" << row << "," << col << ") to file " << options.m_centers << std::endl << std::flush;
    return true;
}

template <typename T>
void Process(MPI_Datatype type) {
	int rank, size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	auto vectorReader = Helper::VectorSetReader::CreateInstance(std::make_shared<Helper::ReaderOptions>(options));
	options.m_inputFiles = Helper::StrUtils::ReplaceAll(options.m_inputFiles, "*", std::to_string(rank));
	if (options.m_inputFiles.find("BIN:") == 0) {
		if (ErrorCode::Success != vectorReader->LoadBinaryFile(options.m_inputFiles.substr(4)))
		{
			fprintf(stderr, "Failed to read input file.\n");
			exit(1);
		}
	}
	else {
		if (ErrorCode::Success != vectorReader->LoadFile(options.m_inputFiles))
		{
			fprintf(stderr, "Failed to read input file.\n");
			exit(1);
		}
	}
    std::shared_ptr<VectorSet> vectors = vectorReader->GetVectorSet();
	std::shared_ptr<MetadataSet> metas = vectorReader->GetMetadataSet();

    COMMON::Dataset<T> data(vectors->Count(), vectors->Dimension(), (T*)vectors->GetData());
    COMMON::KmeansArgs<T> args(options.m_clusterNum, vectors->Dimension(), vectors->Count(), options.m_threadNum);
    std::vector<SizeType> localindices(data.R(), 0);
    for (SizeType i = 0; i < data.R(); i++) localindices[i] = i;

	std::cout << "rank " << rank << " data:(" << data.R() << "," << data.C() << ") machines:" << size << 
		" clusters:" << options.m_clusterNum << " type:" << ((int)options.m_inputValueType) << 
		" threads:" << options.m_threadNum << " lambda:" << options.m_lambda << 
		" samples:" << options.m_localSamples << std::endl << std::flush;
	
    if (rank == 0) {
        std::cout << "rank 0 init centers" << std::endl << std::flush;
		if (!LoadCenters(args.newTCenters, args._K, args._D)) {
			if (options.m_seed >= 0) std::srand(options.m_seed);
			COMMON::InitCenters<T>(data, localindices, 0, data.R(), args, options.m_distMethod, options.m_localSamples, options.m_initIter);
		}
    }

    float currDiff = 1.0, d, currDist, minClusterDist = MaxDist;
    int iteration = 0;
    int noImprovement = 0;
    while (currDiff > options.m_stopDifference && iteration < options.m_maxIter) {
        if (rank == 0) {
            std::memcpy(args.centers, args.newTCenters, sizeof(T)*args._K*args._D);
        }
        MPI_Bcast(args.centers, args._K*args._D, type, 0, MPI_COMM_WORLD);

        args.ClearCenters();
        args.ClearCounts();
        args.ClearDists(-MaxDist);
        d = COMMON::KmeansAssign<T>(data, localindices, 0, data.R(), args, true, options.m_distMethod, (iteration == 0)? 0.0f : options.m_lambda);
		std::memset(args.counts, 0, sizeof(SizeType) * args._K);
		currDist = 0;
		MPI_Allreduce(args.newCounts, args.counts, args._K, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&d, &currDist, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
 
        if (currDist < minClusterDist) {
            noImprovement = 0;
            minClusterDist = currDist;
        }
        else {
            noImprovement++;
        }
        if (noImprovement >= 10) break;

        if (rank == 0) {
            MPI_Reduce(MPI_IN_PLACE, args.newCenters, args._K * args._D, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            currDiff = COMMON::RefineCenters<T>(data, args, options.m_distMethod);
            std::cout << "iter " << iteration << " dist:" << currDist << " diff:" << currDiff << std::endl << std::flush;
        } else
            MPI_Reduce(args.newCenters, args.newCenters, args._K * args._D, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		
        iteration++;
        MPI_Bcast(&currDiff, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
	if (rank == 0) {
		SaveCenters(args.centers, args._K, args._D);
		for (int i = 0; i < args._K; i++)
			std::cout << "cluster " << i << " contains vectors:" << args.counts[i] << std::endl << std::flush;
	}
    MPI_Barrier(MPI_COMM_WORLD);

	if (options.m_outdir.compare("-") != 0) {
		for (int i = 0; i < args._K; i++) {
			if (i % size == rank) {
				std::ofstream out(options.m_outdir + "/" + options.m_outfile + "." + std::to_string(i), std::ios::binary);
				std::ofstream metaout(options.m_outdir + "/" + options.m_outmetafile + "." + std::to_string(i), std::ios::binary);
				std::ofstream metaindexout(options.m_outdir + "/" + options.m_outmetaindexfile + "." + std::to_string(i), std::ios::binary);
				if (!out.is_open() || !metaout.is_open() || !metaindexout.is_open()) {
					std::cout << "Error open write file " << options.m_outfile << " " << options.m_outmetafile << " " << options.m_outmetaindexfile << std::endl << std::flush;
					exit(1);
				}
				out.write((char *)(&args.counts[i]), sizeof(int));
				out.write((char *)(&args._D), sizeof(int));
				if (metas != nullptr) metaindexout.write((char*)(&args.counts[i]), sizeof(int));
				std::uint64_t offset = 0;
				T* recvbuf = args.newTCenters;
				for (int j = 0; j < size; j++) {
					uint64_t offset_before = offset;
					if (j != rank) {
						int recv = 0;
						MPI_Recv(&recv, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						for (int k = 0; k < recv; k++) {
							MPI_Recv(recvbuf, args._D, type, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							out.write((char*)recvbuf, sizeof(T)*args._D);

							if (metas != nullptr) {
								int len;
								MPI_Recv(&len, 1, MPI_INT, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								char *buf = new char[len];
								MPI_Recv(buf, len, MPI_CHAR, j, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								metaout.write(buf, len);
								delete[] buf;
								metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
								offset += len;
							}
						}
						std::cout << "rank " << rank << " <- rank " << j << ":" << recv << " vectors, " << (offset - offset_before) << " bytes meta" << std::endl << std::flush;
					}
					else {
						for (int k = 0; k < data.R(); k++) {
							if (args.label[k] == i) {
								out.write((char*)(data[localindices[k]]), sizeof(T) * args._D);
								if (metas != nullptr) {
									ByteArray meta = metas->GetMetadata(localindices[k]);
									metaout.write((const char*)meta.Data(), meta.Length());
									metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
									offset += meta.Length();
								}
							}
						}
						std::cout << "rank " << rank << " <- rank " << j << ":" << args.newCounts[i] << " vectors, " << (offset - offset_before) << " bytes meta" << std::endl << std::flush;
					}
				}
				if (metas != nullptr) metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
				out.close();
				metaout.close();
				metaindexout.close();
			}
			else {
				int dest = i % size;
				MPI_Send(&args.newCounts[i], 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
				size_t total_len = 0;
				for (int j = 0; j < data.R(); j++) {
					if (args.label[j] == i) {
						MPI_Send(data[localindices[j]], args._D, type, dest, 1, MPI_COMM_WORLD);
						if (metas != nullptr) {
							ByteArray meta = metas->GetMetadata(localindices[j]);
							int len = (int)meta.Length();
							MPI_Send(&len, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
							MPI_Send(meta.Data(), len, MPI_CHAR, dest, 3, MPI_COMM_WORLD);
							total_len += len;
						}
					}
				}
				std::cout << "rank " << rank << " -> rank " << dest << ":" << args.newCounts[i] << " vectors, " << total_len << " bytes meta" << std::endl << std::flush;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
    MPI_Finalize();
}

int main(int argc, char* argv[]) {
    if (!options.Parse(argc - 1, argv + 1))
    {
        exit(1);
    }

    switch (options.m_inputValueType) {
    case SPTAG::VectorValueType::Float:
        Process<float>(MPI_FLOAT);
        break;
    case SPTAG::VectorValueType::Int16:
        Process<std::int16_t>(MPI_SHORT);
        break;
    case SPTAG::VectorValueType::Int8:
        Process<std::int8_t>(MPI_CHAR);
        break;
    case SPTAG::VectorValueType::UInt8:
        Process<std::uint8_t>(MPI_CHAR);
        break;
    default:
        std::cout << "Error data type!" << std::endl;
    }
    return 0;
}
