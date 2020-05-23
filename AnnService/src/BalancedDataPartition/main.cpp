#include <mpi.h>
#include <algorithm>
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
        AddOptionalOption(m_clusterassign, "-a", "--assign", "Number of clusters to be assigned (1<=assign<=4)");
        AddOptionalOption(m_vectorfactor, "-v", "--vectorscale", "Max vector number scale factor.");
        AddOptionalOption(m_closurefactor, "-f", "--closurescale", "Max closure factor");
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
    int m_clusterassign = 1;
    float m_vectorfactor = 2.0f;
    float m_closurefactor = 1.2f;
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
inline float MultipleClustersAssign(const COMMON::Dataset<T>& data,
    std::vector<SizeType>& indices,
    const SizeType first, const SizeType last, COMMON::KmeansArgs<T>& args, bool updateCenters, float lambda) {
    float(*fComputeDistance)(const T* pX, const T* pY, DimensionType length) = COMMON::DistanceCalcSelector<T>(options.m_distMethod);
    float currDist = 0;
    SizeType subsize = (last - first - 1) / args._T + 1;

#pragma omp parallel for num_threads(args._T) shared(data, indices) reduction(+:currDist)
    for (int tid = 0; tid < args._T; tid++)
    {
        SizeType istart = first + tid * subsize;
        SizeType iend = min(first + (tid + 1) * subsize, last);
        SizeType *inewCounts = args.newCounts + tid * args._K;
        float *inewCenters = args.newCenters + tid * args._K * args._D;
        SizeType * iclusterIdx = args.clusterIdx + tid * args._K;
        float * iclusterDist = args.clusterDist + tid * args._K;
        float idist = 0;
        std::vector<SPTAG::COMMON::HeapCell> centerDist(args._K, SPTAG::COMMON::HeapCell());
        for (SizeType i = istart; i < iend; i++) {
            for (int k = 0; k < args._K; k++) {
                float dist = fComputeDistance(data[indices[i]], args.centers + k*args._D, args._D) + lambda*args.counts[k];
                centerDist[k].node = k;
                centerDist[k].distance = dist;
            }
			std::sort(centerDist.begin(), centerDist.end(), [](const COMMON::HeapCell &a, const COMMON::HeapCell &b) {
				return (a.distance < b.distance) || (a.distance == b.distance && a.node < b.node);
			});
            args.label[i] = 0;
            for (int k = 0; k < options.m_clusterassign; k++) {
                if (centerDist[k].distance <= centerDist[0].distance * options.m_closurefactor) {
                    args.label[i] |= ((centerDist[k].node & 0xff) << (k * 8));
                    inewCounts[centerDist[k].node]++;
                    idist += centerDist[k].distance;

                    if (updateCenters) {
                        const T* v = (const T*)data[indices[i]];
                        float* center = inewCenters + centerDist[k].node*args._D;
                        for (DimensionType j = 0; j < args._D; j++) center[j] += v[j];
                        if (centerDist[k].distance > iclusterDist[centerDist[k].node]) {
                            iclusterDist[centerDist[k].node] = centerDist[k].distance;
                            iclusterIdx[centerDist[k].node] = indices[i];
                        }
                    }
                    else {
                        if (centerDist[k].distance <= iclusterDist[centerDist[k].node]) {
                            iclusterDist[centerDist[k].node] = centerDist[k].distance;
                            iclusterIdx[centerDist[k].node] = indices[i];
                        }
                    }
                }
                else {
                    args.label[i] |= (0xff << (k * 8));
                }
            }
            for (int k = options.m_clusterassign; k < 4; k++) {
                args.label[i] |= (0xff << (k * 8));
            }
        }
        currDist += idist;
    }

    for (int i = 1; i < args._T; i++) {
        for (int k = 0; k < args._K; k++)
            args.newCounts[k] += args.newCounts[i*args._K + k];
    }

    if (updateCenters) {
        for (int i = 1; i < args._T; i++) {
            float* currCenter = args.newCenters + i*args._K*args._D;
            for (size_t j = 0; j < ((size_t)args._K) * args._D; j++) args.newCenters[j] += currCenter[j];

            for (int k = 0; k < args._K; k++) {
                if (args.clusterIdx[i*args._K + k] != -1 && args.clusterDist[i*args._K + k] > args.clusterDist[k]) {
                    args.clusterDist[k] = args.clusterDist[i*args._K + k];
                    args.clusterIdx[k] = args.clusterIdx[i*args._K + k];
                }
            }
        }
    }
    else {
        for (int i = 1; i < args._T; i++) {
            for (int k = 0; k < args._K; k++) {
                if (args.clusterIdx[i*args._K + k] != -1 && args.clusterDist[i*args._K + k] <= args.clusterDist[k]) {
                    args.clusterDist[k] = args.clusterDist[i*args._K + k];
                    args.clusterIdx[k] = args.clusterIdx[i*args._K + k];
                }
            }
        }
    }
    return currDist;
}

template <typename T>
inline float HardMultipleClustersAssign(const COMMON::Dataset<T>& data,
	std::vector<SizeType>& indices,
	const SizeType first, const SizeType last, COMMON::KmeansArgs<T>& args, SizeType* mylimit,
	const int clusternum, const bool fill) {
	float(*fComputeDistance)(const T* pX, const T* pY, DimensionType length) = COMMON::DistanceCalcSelector<T>(options.m_distMethod);
	float currDist = 0;
	int threads = 1;
	SizeType subsize = (last - first - 1) / threads + 1;

#pragma omp parallel for num_threads(threads) shared(data, indices) reduction(+:currDist)
	for (int tid = 0; tid < threads; tid++)
	{
		SizeType istart = first + tid * subsize;
		SizeType iend = min(first + (tid + 1) * subsize, last);
		SizeType *inewCounts = args.newCounts + tid * args._K;
		float idist = 0;
		std::vector<SPTAG::COMMON::HeapCell> centerDist(args._K, SPTAG::COMMON::HeapCell());
		for (SizeType i = istart; i < iend; i++) {
			for (int k = 0; k < args._K; k++) {
				float dist = fComputeDistance(data[indices[i]], args.centers + k*args._D, args._D);
				centerDist[k].node = k;
				centerDist[k].distance = dist;
			}
			std::sort(centerDist.begin(), centerDist.end(), [](const COMMON::HeapCell &a, const COMMON::HeapCell &b) {
				return (a.distance < b.distance) || (a.distance == b.distance && a.node < b.node);
			});

			if (centerDist[clusternum].distance <= centerDist[0].distance * options.m_closurefactor &&
				inewCounts[centerDist[clusternum].node] < mylimit[centerDist[clusternum].node]) {
				args.label[i] |= ((centerDist[clusternum].node & 0xff) << (clusternum * 8));
				inewCounts[centerDist[clusternum].node]++;
				idist += centerDist[clusternum].distance;
			}
			else {
				args.label[i] |= (0xff << (clusternum * 8));
			}

			if (fill) {
				for (int k = clusternum + 1; k < 4; k++) {
					args.label[i] |= (0xff << (k * 8));
				}
			}
		}
		currDist += idist;
	}

	std::memset(args.counts, 0, sizeof(SizeType) * args._K);
	for (int i = 0; i < threads; i++) {
		for (int k = 0; k < args._K; k++)
			args.counts[k] += args.newCounts[i*args._K + k];
	}
	return currDist;
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
	unsigned long long localCount = data.R(), totalCount;
	MPI_Allreduce(&localCount, &totalCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	totalCount = static_cast<unsigned long long>(totalCount * 1.0 / args._K * options.m_vectorfactor);

    std::cout << "rank " << rank << " data:(" << data.R() << "," << data.C() << ") machines:" << size << 
        " clusters:" << options.m_clusterNum << " type:" << ((int)options.m_inputValueType) << 
        " threads:" << options.m_threadNum << " lambda:" << options.m_lambda << 
        " samples:" << options.m_localSamples << " maxcountperpartition:" << totalCount << std::endl << std::flush;
    
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
        d = MultipleClustersAssign<T>(data, localindices, 0, data.R(), args, true, (iteration == 0) ? 0.0f : options.m_lambda);
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
    if (options.m_maxIter == 0) {
        if (rank == 0) {
            std::memcpy(args.centers, args.newTCenters, sizeof(T)*args._K*args._D);
        }
        MPI_Bcast(args.centers, args._K*args._D, type, 0, MPI_COMM_WORLD);
    }
    else {
        if (rank == 0) {
            for (int i = 0; i < args._K; i++)
                std::cout << "cluster " << i << " contains vectors:" << args.counts[i] << std::endl << std::flush;
        }
    }
	d = 0;
	std::vector<SizeType> myLimit(args._K, totalCount);
	std::memset(args.counts, 0, sizeof(SizeType)*args._K);
	std::memset(args.label, 0, sizeof(int)*localCount);
    args.ClearCounts();
	for (int i = 0; i < options.m_clusterassign - 1; i++) {
		d += HardMultipleClustersAssign<T>(data, localindices, 0, data.R(), args, myLimit.data(), i, false);
		std::memcpy(myLimit.data(), args.counts, sizeof(SizeType)*args._K);
		MPI_Allreduce(MPI_IN_PLACE, args.counts, args._K, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if (rank == 0) {
			std::cout << "assign " << i << "....................d:" << d << std::endl << std::flush;
			for (int i = 0; i < args._K; i++)
				std::cout << "cluster " << i << " contains vectors:" << args.counts[i] << std::endl << std::flush;
		}
		for (int k = 0; k < args._K; k++)
			myLimit[k] += (totalCount - args.counts[k]) / size;
	}
	d += HardMultipleClustersAssign<T>(data, localindices, 0, data.R(), args, myLimit.data(), options.m_clusterassign - 1, true);
	std::memcpy(args.newCounts, args.counts, sizeof(SizeType)*args._K);

	MPI_Allreduce(args.newCounts, args.counts, args._K, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&d, &currDist, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    if (rank == 0) {
        SaveCenters(args.centers, args._K, args._D);
        std::cout << "final dist:" << currDist << std::endl;
        for (int i = 0; i < args._K; i++)
            std::cout << "cluster " << i << " contains vectors:" << args.counts[i] << std::endl << std::flush;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (options.m_outdir.compare("-") != 0) {
        for (int i = 0; i < args._K; i++) {
            noImprovement = 0;
            if (i % size == rank) {
                std::cout << "Cluster " << i << " start ......" << std::endl << std::flush;
                std::string vecfile = options.m_outdir + "/" + options.m_outfile + "." + std::to_string(i);
                if (fileexists(vecfile.c_str())) noImprovement = 1;
            }
            MPI_Bcast(&noImprovement, 1, MPI_INT, i % size, MPI_COMM_WORLD);
            if (noImprovement) continue;

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
                int recvmetabuflen = 200;
                char* recvmetabuf = new char [recvmetabuflen];
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
                                if (len > recvmetabuflen) {
                                    std::cout << "enlarge recv meta buf to " << len << std::endl << std::flush;
                                    delete[] recvmetabuf;
                                    recvmetabuflen = len;
                                    recvmetabuf = new char[recvmetabuflen];
                                }
                                MPI_Recv(recvmetabuf, len, MPI_CHAR, j, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                metaout.write(recvmetabuf, len);
                                metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
                                offset += len;
                            }
                        }
                        std::cout << "rank " << rank << " <- rank " << j << ":" << recv << " vectors, " << (offset - offset_before) << " bytes meta" << std::endl << std::flush;
                    }
                    else {
                        size_t total_rec = 0;
                        for (int k = 0; k < data.R(); k++) {
                            for (int kk = 0; kk < 32; kk += 8) {
                                if (((args.label[k] >> kk) & 0xff) == (i & 0xff)) {
                                    out.write((char*)(data[localindices[k]]), sizeof(T) * args._D);
                                    if (metas != nullptr) {
                                        ByteArray meta = metas->GetMetadata(localindices[k]);
                                        metaout.write((const char*)meta.Data(), meta.Length());
                                        metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
                                        offset += meta.Length();
                                    }
                                    total_rec++;
                                }
                            }
                        }
                        std::cout << "rank " << rank << " <- rank " << j << ":" << args.newCounts[i] << "(" << total_rec << ") vectors, " << (offset - offset_before) << " bytes meta" << std::endl << std::flush;
                    }
                }
                delete[] recvmetabuf;
                if (metas != nullptr) metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
                out.close();
                metaout.close();
                metaindexout.close();
            }
            else {
                int dest = i % size;
                MPI_Send(&args.newCounts[i], 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                size_t total_len = 0;
                size_t total_rec = 0;
                for (int j = 0; j < data.R(); j++) {
                    for (int kk = 0; kk < 32; kk += 8) {
                        if (((args.label[j] >> kk) & 0xff) == (i & 0xff)) {
                            MPI_Send(data[localindices[j]], args._D, type, dest, 1, MPI_COMM_WORLD);
                            if (metas != nullptr) {
                                ByteArray meta = metas->GetMetadata(localindices[j]);
                                int len = (int)meta.Length();
                                MPI_Send(&len, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
                                MPI_Send(meta.Data(), len, MPI_CHAR, dest, 3, MPI_COMM_WORLD);
                                total_len += len;
                            }
                            total_rec++;
                        }
                    }
                }
                std::cout << "rank " << rank << " -> rank " << dest << ":" << args.newCounts[i] << "(" << total_rec << ") vectors, " << total_len << " bytes meta" << std::endl << std::flush;
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
