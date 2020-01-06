#include <mpi.h>
#include "inc/Core/Common/DistanceUtils.h"
#include "inc/Core/Common/Dataset.h"
#include "inc/Core/Common/BKTree.h"
#include "inc/Helper/VectorSetReader.h"

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
	}

	~PartitionOptions() {}

	std::string m_inputFiles;
	int m_clusterNum;
	
	float m_stopDifference = 0.0001f;
	int m_maxIter = 100;
	int m_localSamples = 1000;
	float m_lambda = 0.000001f;
	DistCalcMethod m_distMethod = DistCalcMethod::L2;
	
	std::string m_outfile = "vectors.bin";
	std::string m_outmetafile = "meta.bin";
	std::string m_outmetaindexfile = "metaindex.bin";
};

PartitionOptions options;

template <typename T>
void Process(MPI_Datatype type) {
	auto vectorReader = Helper::VectorSetReader::CreateInstance(std::make_shared<Helper::ReaderOptions>(options));
	if (ErrorCode::Success != vectorReader->LoadFile(options.m_inputFiles))
	{
		fprintf(stderr, "Failed to read input file.\n");
		exit(1);
	}

	std::shared_ptr<VectorSet> vectors = vectorReader->GetVectorSet();
	std::shared_ptr<MetadataSet> metas = vectorReader->GetMetadataSet();

	COMMON::Dataset<T> data(vectors->Count(), vectors->Dimension(), (T*)vectors->GetData());
	std::vector<SizeType> localindices;
	localindices.resize(data.R());
	for (SizeType i = 0; i < localindices.size(); i++) localindices[i] = i;

	int rank, size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::cout << "data:(" << data.R() << "," << data.C() << ") machines:" << size << " clusters:" << options.m_clusterNum << std::endl;

        std::cout << "clusterNum:" << options.m_clusterNum << " dimension:" << vectors->Dimension() << " type:" << ((int)options.m_inputValueType) << " threads:" << options.m_threadNum << std::endl << std::flush;
	COMMON::KmeansArgs<T> args(options.m_clusterNum, vectors->Dimension(), vectors->Count(), options.m_threadNum);

	if (rank == 0) {
		std::cout << "rank 0 init centers" << std::endl << std::flush;
		COMMON::InitCenters<T>(data, localindices, 0, data.R(), args, options.m_distMethod, options.m_localSamples);
	}

	float currDiff = 1.0, currDist, minClusterDist = MaxDist;
        int iteration = 0;
	int noImprovement = 0;
	while (currDiff > options.m_stopDifference && iteration < options.m_maxIter) {
		std::memcpy(args.centers, args.newTCenters, sizeof(T)*args._K*data.C());
		MPI_Bcast(args.centers, args._K*data.C(), type, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			std::cout << "iter " << iteration << " rank 0 bcast centers" << std::endl << std::flush;
		}
		
		args.ClearCenters();
		args.ClearCounts();
		args.ClearDists(-MaxDist);
		float d = COMMON::KmeansAssign<T>(data, localindices, 0, data.R(), args, true, options.m_distMethod);
		MPI_Allreduce(args.newCounts, args.counts, args._K, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&d, &currDist, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
 
		if (currDist < minClusterDist) {
			noImprovement = 0;
			minClusterDist = currDist;
		}
		else {
			noImprovement++;
		}
		if (noImprovement >= 5) break;

                if (rank == 0) {
                    MPI_Reduce(MPI_IN_PLACE, args.newCenters, args._K * data.C(), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	            currDiff = COMMON::RefineCenters<T>(data, args, options.m_distMethod);
                } else
                    MPI_Reduce(args.newCenters, args.newCenters, args._K * data.C(), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
 
                iteration++;
		MPI_Bcast(&currDiff, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}

	std::cout << "begin write data" << std::endl << std::flush;
	MPI_Request* req = new MPI_Request[args._K];
	MPI_Status* sta = new MPI_Status[args._K];

	for (int i = 0; i < args._K; i++) {
		if (i == rank) continue;

		T* send = new T[args.newCounts[i] * data.C()];
		int k = 0;
		int len = 0;
		for (int j = 0; j < data.R(); j++) {
			if (args.label[j] == i) {
				std::memcpy(send + k * data.C(), data[j], sizeof(T) * data.C());
				k++;
				len += (int)(metas->GetMetadata(j).Length()) + 1;
			}
		}
		MPI_Isend(&args.newCounts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req[rank]);
		MPI_Wait(&req[rank], &sta[rank]);
		MPI_Isend(send, args.newCounts[i] * data.C(), type, i, 1, MPI_COMM_WORLD, &req[rank]);
		MPI_Wait(&req[rank], &sta[rank]);
		delete[] send;

		char* buf = new char[len];
		k = 0;
		for (int j = 0; j < data.R(); j++) {
			if (args.label[j] == i) {
				ByteArray meta = metas->GetMetadata(j);
				std::memcpy(buf + k, meta.Data(), meta.Length());
				k += (int)meta.Length();
				buf[k++] = '\n';
			}
		}
		MPI_Isend(&k, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &req[rank]);
		MPI_Wait(&req[rank], &sta[rank]);
		MPI_Isend(buf, k, MPI_CHAR, i, 3, MPI_COMM_WORLD, &req[rank]);
		MPI_Wait(&req[rank], &sta[rank]);
		delete[] buf;
	}

	std::cout << "rank " << rank << " begin to recv data" << std::endl;

	std::ofstream out(options.m_outfile, std::ios::binary);
	std::ofstream metaout(options.m_outmetafile, std::ios::binary);
	std::ofstream metaindexout(options.m_outmetaindexfile, std::ios::binary);
	if (!out.is_open() || !metaout.is_open() || !metaindexout.is_open()) {
		std::cout << "Error open write file " << options.m_outfile << " " << options.m_outmetafile << " " << options.m_outmetaindexfile << std::endl;
	}
	else {
		out.write((char *)(&args.counts[rank]), sizeof(int));
		out.write((char *)(&args._D), sizeof(int));
		metaindexout.write((char*)(&args.counts[rank]), sizeof(int));
		std::uint64_t offset = 0;
		for (int i = 0; i < args._K; i++) {
			if (i != rank) {
				int recv;
				MPI_Irecv(&recv, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req[i]);
				MPI_Wait(&req[i], &sta[i]);
				T *recvbuf = new T[recv * args._D];
				MPI_Irecv(recvbuf, recv * args._D, type, i, 1, MPI_COMM_WORLD, &req[i]);
				MPI_Wait(&req[i], &sta[i]);
				out.write((char*)recvbuf, sizeof(T)*recv*args._D);
				delete[] recvbuf;

				MPI_Irecv(&recv, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &req[i]);
				MPI_Wait(&req[i], &sta[i]);
				char *buf = new char[recv];
				MPI_Irecv(buf, recv, MPI_CHAR, i, 3, MPI_COMM_WORLD, &req[i]);
				MPI_Wait(&req[i], &sta[i]);

				char* current;
				char* context = NULL;
				current = strtok_s(buf, "\n", &context);
				while (current != NULL) {
					int slen = (int)(strlen(current));
					metaout.write(current, slen);
					metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
					offset += slen;
					current = strtok_s(NULL, "\n", &context);
				}
				delete[] buf;

			}
			else {
				for (int j = 0; j < data.R(); j++) {
					if (args.label[j] == i) {
						out.write((char*)(data[j]), sizeof(T) * args._D);
						ByteArray meta = metas->GetMetadata(j);
						metaout.write((const char*)meta.Data(), meta.Length());
						metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
						offset += meta.Length();
					}
				}
			}
		}
		out.close();
		metaout.close();
		metaindexout.close();
	}
	delete[] req;
	delete[] sta;
	MPI_Finalize();
}

int main(int argc, char* argv[]) {
	if (!options.Parse(argc - 1, argv + 1))
	{
		exit(1);
	}

	omp_set_num_threads(options.m_threadNum);

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
