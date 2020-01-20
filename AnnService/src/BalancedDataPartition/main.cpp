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
} options;

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
    COMMON::KmeansArgs<T> args(options.m_clusterNum, vectors->Dimension(), vectors->Count(), options.m_threadNum);
    std::vector<SizeType> localindices(data.R(), 0);
    for (SizeType i = 0; i < data.R(); i++) localindices[i] = i;

    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "rank " << rank << " data:(" << data.R() << "," << data.C() << ") machines:" << size << " clusters:" << options.m_clusterNum << " type:" << ((int)options.m_inputValueType) << " threads:" << options.m_threadNum << " lambda:" << options.m_lambda << std::endl << std::flush;

    if (rank == 0) {
        std::cout << "rank 0 init centers" << std::endl << std::flush;
        COMMON::InitCenters<T>(data, localindices, 0, data.R(), args, options.m_distMethod, options.m_localSamples);
    }

    float currDiff = 1.0, currDist, minClusterDist = MaxDist;
    int iteration = 0;
    int noImprovement = 0;
    while (currDiff > options.m_stopDifference && iteration < options.m_maxIter) {
        if (rank == 0) {
            std::cout << "iter " << iteration << " rank 0 bcast centers" << std::endl << std::flush;
            std::memcpy(args.centers, args.newTCenters, sizeof(T)*args._K*args._D);
        }
        MPI_Bcast(args.centers, args._K*args._D, type, 0, MPI_COMM_WORLD);

        args.ClearCenters();
        args.ClearCounts();
        args.ClearDists(-MaxDist);
        float d = COMMON::KmeansAssign<T>(data, localindices, 0, data.R(), args, true, options.m_distMethod, options.m_lambda);
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
            MPI_Reduce(MPI_IN_PLACE, args.newCenters, args._K * args._D, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            currDiff = COMMON::RefineCenters<T>(data, args, options.m_distMethod);
        } else
            MPI_Reduce(args.newCenters, args.newCenters, args._K * args._D, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
 
        iteration++;
        MPI_Bcast(&currDiff, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    std::cout << "rank " << rank << " contains vectors:" << args.counts[rank] << std::endl << std::flush;

    std::cout << "begin write data" << std::endl << std::flush;
    std::ofstream out(options.m_outfile, std::ios::binary);
    std::ofstream metaout(options.m_outmetafile, std::ios::binary);
    std::ofstream metaindexout(options.m_outmetaindexfile, std::ios::binary);
    if (!out.is_open() || !metaout.is_open() || !metaindexout.is_open()) {
        std::cout << "Error open write file " << options.m_outfile << " " << options.m_outmetafile << " " << options.m_outmetaindexfile << std::endl << std::flush;
        exit(1);
    }
    for (int i = 0; i < args._K; i++) {
        if (i == rank) {
            out.write((char *)(&args.counts[rank]), sizeof(int));
            out.write((char *)(&args._D), sizeof(int));
            metaindexout.write((char*)(&args.counts[rank]), sizeof(int));
            std::uint64_t offset = 0;
            T* recvbuf = args.newTCenters;
            for (int j = 0; j < args._K; j++) {
                if (j != rank) {
                    int recv = 0;
                    size_t total_len = 0;
                    MPI_Recv(&recv, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (int k = 0; k < recv; k++) {
                        MPI_Recv(recvbuf, args._D, type, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        out.write((char*)recvbuf, sizeof(T)*args._D);

                        int len;
                        MPI_Recv(&len, 1, MPI_INT, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        char *buf = new char[len];
                        MPI_Recv(buf, len, MPI_CHAR, j, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        metaout.write(buf, len);
                        delete[] buf;
                        metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
                        offset += len;
                        total_len += len;
                    }
                    std::cout << "rank " << rank << " <- rank " << j << ":" << recv << " vectors, " << total_len << " bytes meta" << std::endl << std::flush;
                }
                else {
                    for (int k = 0; k < data.R(); k++) {
                        if (args.label[k] == j) {
                            out.write((char*)(data[localindices[k]]), sizeof(T) * args._D);
                            ByteArray meta = metas->GetMetadata(localindices[k]);
                            metaout.write((const char*)meta.Data(), meta.Length());
                            metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
                            offset += meta.Length();
                        }
                    }
                }
            }
            metaindexout.write((char*)(&offset), sizeof(std::uint64_t));
        } else {
            MPI_Send(&args.newCounts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            size_t total_len = 0;
            for (int j = 0; j < data.R(); j++) {
                if (args.label[j] == i) {
                    MPI_Send(data[localindices[j]], args._D, type, i, 1, MPI_COMM_WORLD);
                    ByteArray meta = metas->GetMetadata(localindices[j]);
                    int len = meta.Length();
                    MPI_Send(&len, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                    MPI_Send(meta.Data(), len, MPI_CHAR, i, 3, MPI_COMM_WORLD);
                    total_len += len;
                }
            }
            std::cout << "rank " << rank << " -> rank " << i << ":" << args.newCounts[i] << " vectors, " << total_len << " bytes meta" << std::endl << std::flush;
        }
    }
    out.close();
    metaout.close();
    metaindexout.close();
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
