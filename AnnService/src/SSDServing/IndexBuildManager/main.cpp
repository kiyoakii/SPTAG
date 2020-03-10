#include "inc/SSDServing/Common/stdafx.h"
#include "inc/SSDServing/SelectHead_BKT/BootSelectHead.h"
#include "inc/SSDServing/SelectHead_BKT/Options.h"
#include "inc/SSDServing/BuildHead/BootBuildHead.h"
#include "inc/SSDServing/BuildHead/Options.h"
#include "inc/SSDServing/VectorSearch/BootVectorSearch.h"
#include "inc/SSDServing/VectorSearch/Options.h"
using namespace SPTAG;

namespace SPTAG {
	namespace SSDServing {
		int internalMain(int argc, char* argv[]) {
			if (argc < 2)
			{
				fprintf(stderr,
					"SSDServing.exe configFilePath\n");
			}

			Helper::IniReader iniReader;
			iniReader.LoadIniFile(argv[1]);

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

			return 0;
		}
	}
}