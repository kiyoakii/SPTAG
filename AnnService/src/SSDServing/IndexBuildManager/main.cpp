#include "inc/Helper/SimpleIniReader.h"

#include "inc/SSDServing/IndexBuildManager/CommonDefines.h"
#include "inc/SSDServing/IndexBuildManager/Options.h"
#include "inc/SSDServing/SelectHead_BKT/BootSelectHead.h"
#include "inc/SSDServing/SelectHead_BKT/Options.h"
#include "inc/SSDServing/BuildHead/BootBuildHead.h"
#include "inc/SSDServing/BuildHead/Options.h"
#include "inc/SSDServing/VectorSearch/BootVectorSearch.h"
#include "inc/SSDServing/VectorSearch/Options.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"

using namespace SPTAG;

namespace SPTAG {
	namespace SSDServing {
		
		BaseOptions COMMON_OPTS;

		int internalMain(int argc, char* argv[]) {
			if (argc < 2)
			{
				fprintf(stderr,
					"ssdserving configFilePath\n");
				exit(-1);
			}

			Helper::IniReader iniReader;
			iniReader.LoadIniFile(argv[1]);

			VectorSearch::TimeUtils::StopW sw;

			auto& baseParameters = iniReader.GetParameters("Base");
			if (!baseParameters.empty())
			{
				for (const auto& iter : baseParameters)
				{
					COMMON_OPTS.SetParameter(iter.first.c_str(), iter.second.c_str());
				}
			}

			auto& selectHeadParameters = iniReader.GetParameters("SelectHead");
			if (!selectHeadParameters.empty())
			{
				SSDServing::SelectHead_BKT::Options slOpts;
				for (const auto& iter : selectHeadParameters)
				{
					slOpts.SetParameter(iter.first.c_str(), iter.second.c_str());
				}
				SSDServing::SelectHead_BKT::Bootstrap(slOpts);
			}
			double selectHeadTime = sw.getElapsedSec();
			sw.reset();

			auto& buildHeadParameters = iniReader.GetParameters("BuildHead");
			if (!buildHeadParameters.empty())
			{
				SSDServing::BuildHead::Options bhOpts;
				for (const auto& iter : buildHeadParameters)
				{
					bhOpts.SetParameter(iter.first.c_str(), iter.second.c_str());
				}
				SSDServing::BuildHead::Bootstrap(bhOpts);
			}
			double buildHeadTime = sw.getElapsedSec();
			sw.reset();

			auto& buildSSDParameters = iniReader.GetParameters("BuildSSDIndex");
			if (!buildSSDParameters.empty())
			{
				SSDServing::VectorSearch::Options vsOpts;
				for (const auto& iter : buildSSDParameters)
				{
					vsOpts.SetParameter(iter.first.c_str(), iter.second.c_str());
				}
				SSDServing::VectorSearch::Bootstrap(vsOpts);
			}
			double buildSSDTime = sw.getElapsedSec();
			sw.reset();

			auto& searchSSDParameters = iniReader.GetParameters("SearchSSDIndex");
			if (!searchSSDParameters.empty())
			{
				SSDServing::VectorSearch::Options vsOpts;
				for (const auto& iter : searchSSDParameters)
				{
					vsOpts.SetParameter(iter.first.c_str(), iter.second.c_str());
				}
				SSDServing::VectorSearch::Bootstrap(vsOpts);
			}
			double searchSSDTime = sw.getElapsedSec();

			fprintf(stderr, "select head time: %.2lf\nbuild head time: %.2lf\nbuild ssd time: %.2lf\nsearch ssd time: %.2lf\n", 
				selectHeadTime,
				buildHeadTime,
				buildSSDTime,
				searchSSDTime
			);

			exit(0);
		}
	}
}

// switch between exe and static library by _$(OutputType)
#ifdef _exe

int main(int argc, char* argv[]) {
	SPTAG::SSDServing::internalMain(argc, argv);
}

#endif