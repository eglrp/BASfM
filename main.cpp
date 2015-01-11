#include "BAutil.h"


void Help()
{
	printf("Usage: BAfromVSFM.exe <Path> <Fixed intrinsic NVM> <nvm file> <Shared intrinsic cameraInfo> <cameraInfo file> <optim shared intrinisc>  \n <Convert Lens model> <Inlier threshold><Intrinsic file>\n");
	printf("<nvm file>: result from Visual SFM\n");
	printf("<camera ini file>: filename width height  intrinisc options, lens distortion options, extrinsic options, availbiltiy (not needed for shared intrinisc)\n");
	printf("<intrinsic param file>: (optional)  filename fx fy skew u0 v0 radiral(1,2,3) tangential(1,2) prism(1,2)\n");

	printf("Options for intrinisc:\n BA_OPT_INTRINSIC_ALL = 0\n BA_OPT_INTRINSIC_ALL_FIXED = 1\n BA_OPT_INTRINSIC_SKEW_FIXED = 2\n BA_OPT_INTRINSIC_SKEW_ZERO_FIXED = 3\n BA_OPT_INTRINSIC_CENTER_FIXED = 4\n");
	printf("Options for lens distorion:\n BA_LENS_FISHEYE = -1\n BA_OPT_LENSDIST_RADIAL_AND_TANGENTIAL  = 0\n BA_OPT_LENSDIST_ALL_FIXED = 1\n BA_OPT_LENSDIST_RADIAL_1ST_ONLY = 2\n BA_OPT_LENSDIST_RADIAL_ONLY =3\n");
	printf("Options for extrinisc:\n BA_OPT_EXTRINSIC_ALL  = 0\n BA_OPT_EXTRINSIC_ALL_FIXED = 1\n BA_OPT_EXTRINSIC_R_FIXED = 2\n BA_OPT_EXTRINSIC_T_FIXED = 3\n");

	return;
}
int main(int argc, char* argv[])
{

	/*if (argc == 1)
	{
		Help();
		return 0;
	}*/

	argv[1] = "C:/Data/CellPhone/Corpus/VisSfm",
		argv[2] = "0", // fixed intrinsic in NVM
		argv[3] = "Corpus.nvm",
		argv[4] = "0", //Shared cameraInfo
		argv[5] = "CamInfo.txt",
		argv[6] = "2", //Optimize for shared intrinisc
		argv[7] = "0",//Lens conversion
		argv[8] = "5";//thresholdd. If minus value, all points from NVM are regarded as inliers.
		//argv[9] = "Intrinsics.txt",
		argc = 9;

	int NVMfixIntrinisc = atoi(argv[2]),
		sharedCamInfo = atoi(argv[4]),
		sharedIntrinsics = atoi(argv[6]),
		LensConversion = atoi(argv[7]);
	double thresh = atof(argv[8]);

	string rootfolder(argv[1]);

	//---------load an NVM file--------------
	cout << "Loading the NVM file...";
	string nvmfile = argv[3];
	string nvmpath = rootfolder + "/" + nvmfile;

	BA::NVM    nvmdata;
	if (!BA::loadNVM(nvmpath, nvmdata, NVMfixIntrinisc))
		return 1;
	cout << "Done." << endl;
	//-------------------------------------

	//---------convert NVM data to CameraData---------
	cout << "Converting NVM data to CameraData...";
	string imginfofile = rootfolder + "/" + argv[5];
	vector<BA::CameraData> camera_before;
	if (!BA::initCameraData(nvmdata, imginfofile, camera_before, sharedCamInfo))
		return 1;
	cout << "Done." << endl;
	//-------------------------------------


	//---------load initial parameters if applicable---------
	bool loadIntrinsicSucces = false;
	if (argc == 10)
	{
		cout << "Loading initial intrinsic parameters...";
		string intrinsicfile = rootfolder + "/" + argv[9];
		if (!BA::loadInitialIntrinsics(intrinsicfile, nvmdata.filename_id, camera_before, sharedIntrinsics))
			printf("Cannot load %s\n", intrinsicfile);
		else
		{
			cout << "Done." << endl;
			loadIntrinsicSucces = true;
		}
	}
	//-------------------------------------

	//---------perform bundle adjustment---------
	vector< vector<double> > xyz_before(nvmdata.xyz);
	vector< vector<double> > xyz_after(nvmdata.xyz);
	vector<BA::CameraData> camera_after(camera_before);

	if (LensConversion == 1 && loadIntrinsicSucces == false)
	{
		cout << "\n" << "Converting Lens Model...";
		VisSfMLens2OPENCV(camera_before, xyz_before);
		cout << "Done." << endl;
	}

	//---------reprojection error before BA---------
	cout << "Calculating reprojection erros (before) and Saving the result...\n";
	BA::residualData resdata_before;
	BA::calcReprojectionError(camera_before, xyz_before, resdata_before, thresh);
	BA::saveAllData(rootfolder, camera_before, xyz_before, resdata_before, "BA_", false);

	printf("Before (%d/%d):\n Max & Min Reprojection Errors: (%.2f, %.2f) (%.2f, %.2f)\n Mean Absolute Reprojection Errors :  (%.3f %.3f)\n STD Absolute Reprojection Errors : (%.3f %.3f)\n",
		resdata_before.inliersCount, resdata_before.MaxCount,
		resdata_before.MaxXY[0], resdata_before.MaxXY[1], resdata_before.MinXY[0], resdata_before.MinXY[1],
		resdata_before.MeanABSXY[0], resdata_before.MeanABSXY[1],
		resdata_before.STDXY[0], resdata_before.STDXY[1]);

	ceres::Solver::Summary summary;
	ceres::Solver::Options options;
	BA::setCeresOption(nvmdata, options);
	//options.use_nonmonotonic_steps = false;


	cout << "\n" << "Run Bundle Adjustment..." << endl;
	cout << "\t# of Cameras: " << camera_before.size() << "\n" << "\t# of Points : " << nvmdata.n3dPoints << "\n" << endl;

	if (sharedIntrinsics==1)
	{
		BA::runBundleAdjustmentSharedIntrinsic(camera_after, xyz_after, options, summary, thresh);
		for (int ii = 1; ii < camera_after.size(); ii++)
			copyIntrinsic(camera_after[0], camera_after[ii]);
	}
	else if (sharedIntrinsics == 2)
	{
		BA::runBundleAdjustmentPartiallySharedIntrinsic(camera_after, xyz_after, options, summary, thresh);

		BA::CameraData MasterSharedIntrinisc;
		for (int ii = 0; ii < camera_after.size(); ii++)
		{
			if (camera_after[ii].available && camera_after[ii].sharedIntrinisc)
			{
				printf("found the master at %d\n", ii);
				copyIntrinsic(camera_after[ii], MasterSharedIntrinisc);
				break;
			}
		}
	
		for (int ii = 0; ii < camera_after.size(); ii++)
			if (camera_after[ii].available && camera_after[ii].sharedIntrinisc)
				copyIntrinsic(MasterSharedIntrinisc, camera_after[ii]);
	}
	else
		BA::runBundleAdjustment(camera_after, xyz_after, options, summary, thresh);

	cout << summary.BriefReport() << "\n";

	string ceres_report = rootfolder + "/BA_CeresReport.txt";
	ofstream ofs(ceres_report);
	if (ofs.fail())
		cerr << "Cannot write " << ceres_report << endl;
	else
		ofs << summary.FullReport();
	ofs.close();
	//-------------------------------------


	//--------save as NVM format------------------
	cout << "\n";
	cout << "Saving the result as NVM format...";
	BA::saveNVM(rootfolder, nvmfile, camera_after, xyz_after, nvmdata);
	cout << "Done." << endl;
	//-------------------------------------

	//---------reprojection error after BA---------
	cout << "Calculating reprojection erros ( after) and Saving the result...";
	//Caculate the error on inliers before set.
	for (int ii = 0; ii < nvmdata.nCamera; ii++)
		for (int jj = 0; jj < camera_after[ii].inlier.size(); jj++)
			camera_after[ii].inlier[jj] = camera_before[ii].inlier[jj];

	BA::residualData resdata_after;
	BA::calcReprojectionError(camera_after, xyz_after, resdata_after, thresh);
	BA::saveAllData(rootfolder, camera_after, xyz_after, resdata_after, "BA_", true);
	cout << "Done." << endl;
	//-------------------------------------

	printf("After (%d,%d):\n Max & Min Reprojection Errors: (%.2f, %.2f) (%.2f, %.2f)\n Mean Absolute Reprojection Errors : (%.3f %.3f)\n STD Absolute Reprojection Errors : (%.3f %.3f)\n",
		resdata_after.inliersCount, resdata_after.MaxCount,
		resdata_after.MaxXY[0], resdata_after.MaxXY[1], resdata_after.MinXY[0], resdata_after.MinXY[1],
		resdata_after.MeanABSXY[0], resdata_after.MeanABSXY[1],
		resdata_after.STDXY[0], resdata_after.STDXY[1]);

	return 0;
}
