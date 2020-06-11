#include <inc/SSDServing/IndexBuildManager/Utils.h>

SPTAG::SSDServing::Neighbor::Neighbor(SPTAG::SizeType k, float d) : key(k), dist(d) {}

SPTAG::SSDServing::Neighbor::Neighbor(const SPTAG::SSDServing::Neighbor& o) : key(o.key), dist(o.dist) {}

bool SPTAG::SSDServing::Neighbor::operator < (const SPTAG::SSDServing::Neighbor& another) const
{
	return this->dist == another.dist ? this->key < another.key : this->dist < another.dist;
}

void SPTAG::SSDServing::writeTruthFile(const std::string truthFile, size_t queryNumber, const int K, std::vector<std::vector<SPTAG::SizeType>>& truthset, SPTAG::TruthFileType TFT) {

	if (TFT == SPTAG::TruthFileType::TXT)
	{
		std::ofstream of(truthFile);
		for (size_t i = 0; i < queryNumber; i++)
		{
			for (size_t k = 0; k < K; k++)
			{
				of << truthset[i][k];
				if (k != K - 1)
				{
					of << " ";
				}
			}
			of << std::endl;
		}
	}
	else if (TFT == SPTAG::TruthFileType::XVEC)
	{
		std::ofstream of(truthFile, std::ios_base::binary);
		for (size_t i = 0; i < queryNumber; i++)
		{
			of.write(reinterpret_cast<const char*>(&K), 4);
			of.write(reinterpret_cast<char*>(truthset[i].data()), K * 4);
		}
	}
	else {
		fprintf(stderr, "Found unsupported file type for generating truth. ");
		exit(-1);
	}
}
