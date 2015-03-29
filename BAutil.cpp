#include "BAutil.h"

namespace BA
{
	template <typename T>	void PinholeReprojection(const T* const FocalLength, const T* const OpticalCenter, const T* const Skew, const T* const Radialfirst, const T* const Radialothers,
		const T* const Tangential, const T* const Prism, const T* const AngleAxis, const T* const Translation, const T* const point3D, T* reprojected2D)
	{

		T fx = FocalLength[0];
		T fy = FocalLength[1];
		T s = Skew[0];
		T u0 = OpticalCenter[0];
		T v0 = OpticalCenter[1];
		T a0 = Radialfirst[0], a1 = Radialothers[0], a2 = Radialothers[1];
		T p0 = Tangential[0], p1 = Tangential[1];
		T s0 = Prism[0], s1 = Prism[1];

		//rotate 3d point
		T pt[3];
		ceres::AngleAxisRotatePoint(AngleAxis, point3D, pt);

		//translate 3d point
		pt[0] += Translation[0];
		pt[1] += Translation[1];
		pt[2] += Translation[2];

		// Project to normalize coordinate
		T xn = pt[0] / pt[2];
		T yn = pt[1] / pt[2];

		//apply lens distortion
		T r2 = xn*xn + yn*yn;
		T r4 = r2*r2;
		T r6 = r2*r4;
		T X2 = xn*xn;
		T Y2 = yn*yn;
		T XY = xn*yn;

		T radial = T(1.0) + a0*r2 + a1*r4 + a2*r6;
		T tangential_x = p0*(r2 + T(2.0)*X2) + T(2.0)*p1*XY;
		T tangential_y = p1*(r2 + T(2.0)*Y2) + T(2.0)*p0*XY;
		T prism_x = s0*r2;
		T prism_y = s1*r2;

		T xm = radial*xn + tangential_x + prism_x;
		T ym = radial*yn + tangential_y + prism_y;

		reprojected2D[0] = fx*xm + s*ym + u0;
		reprojected2D[1] = fy*ym + v0;
	}
	template <typename T>	void FOVReprojection(const T* const FocalLength, const T* const OpticalCenter, const T* const Skew, const T* const omega, const T* const DistCtr,
		const T* const AngleAxis, const T* const Translation, const T* const point3D, T* reprojected2D)
	{
		T fx = FocalLength[0], fy = FocalLength[1], s = Skew[0], u0 = OpticalCenter[0], v0 = OpticalCenter[1];

		//rotate 3d point
		T pt[3];
		ceres::AngleAxisRotatePoint(AngleAxis, point3D, pt);

		//translate 3d point
		pt[0] += Translation[0], pt[1] += Translation[1], pt[2] += Translation[2];

		// Project to normalize coordinate
		T xn = pt[0] / pt[2], yn = pt[1] / pt[2];

		// Convert to pixel image coord
		T u = fx*xn + s*yn + u0, v = fy*yn + v0;

		//apply lens distortion
		T x = u - DistCtr[0], y = v - DistCtr[1];
		T ru = sqrt(x*x + y*y), rd = atan(2.0*ru*tan(T(0.5)*omega[0])) / omega[0];
		T t = rd / ru;
		T x_u = t*x, y_u = t*y;

		reprojected2D[0] = x_u + DistCtr[0];
		reprojected2D[1] = y_u + DistCtr[1];
	}

	void calcReprojectionError(const vector<CameraData> &camera, const vector< vector<double> > &xyz, residualData &res, double thresh)
	{
		res.error.clear(); res.error.reserve(xyz.size() * 50);
		res.ID.clear(); res.ID.reserve(xyz.size() * 50);
		res.reprojected_pt.reserve(xyz.size() * 50);
		res.observed_pt.reserve(xyz.size() * 50);

		double MeanABSXY[2] = { 0.0, 0.0 }, max_abs_error[2] = { 0.0, 0.0 };
		int count = 0;

		res.MaxCount = 0;
		const int nCameras = camera.size();
		for (int camID = 0; camID < nCameras; camID++)
		{
			const CameraData *cam = &camera[camID];
			if (!cam->available)
				continue;

			const double *focal = cam->FocalLength, *imcenter = cam->OpticalCenter, *skew = &cam->Skew, *rad_1st = &cam->Radialfirst,
				*rad_other = cam->Radialothers, *tangential = cam->Tangential, *prism = cam->Prism,
				*omega = &cam->omega, *DistCtr = cam->DistCtr,
				*angleaxis = cam->AngleAxis, *trans = cam->Translation;

			const int nPoints = cam->point2D.size();
			for (int i = 0; i < nPoints; i++)
			{
				if (!cam->inlier[i])
					continue;

				double observed_x = cam->point2D[i][0], observed_y = cam->point2D[i][1];

				int ptID = cam->ptID[i];
				const double *pt3D = &xyz[ptID][0];

				double reprojected2D[2];
				if (cam->lenstype == BA_LENS_PINHOLE)
					PinholeReprojection(focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D, reprojected2D);
				else
					FOVReprojection(focal, imcenter, skew, omega, DistCtr, angleaxis, trans, pt3D, reprojected2D);

				vector<double> residual(2);
				residual[0] = reprojected2D[0] - observed_x, residual[1] = reprojected2D[1] - observed_y;

				res.MaxCount++;
				if (thresh > 0.0 && sqrt(pow(residual[0], 2) + pow(residual[1], 2)) > thresh)
					continue;

				res.error.push_back(residual);

				if (abs(residual[0]) > max_abs_error[0])
					max_abs_error[0] = abs(residual[0]);
				if (abs(residual[1]) > max_abs_error[1])
					max_abs_error[1] = abs(residual[1]);

				MeanABSXY[0] += abs(residual[0]), MeanABSXY[1] += abs(residual[1]);
				count++;

				vector<int> ID(2);
				ID[0] = ptID;
				ID[1] = camID;
				res.ID.push_back(ID);

				vector<double> reprojected_pt(2);
				reprojected_pt[0] = reprojected2D[0], reprojected_pt[1] = reprojected2D[1];
				res.reprojected_pt.push_back(reprojected_pt);

				res.observed_pt.push_back(cam->point2D[i]);
			}
		}

		res.inliersCount = count;
		vector<double> ReProjectionErrorX, ReProjectionErrorY;
		ReProjectionErrorX.reserve(count), ReProjectionErrorY.reserve(count);
		for (int ii = 0; ii < res.reprojected_pt.size(); ii++)
		{
			ReProjectionErrorX.push_back(abs(res.error[ii][0]));
			ReProjectionErrorY.push_back(abs(res.error[ii][1]));
		}

		res.MinXY[0] = *min_element(ReProjectionErrorX.begin(), ReProjectionErrorX.end());
		res.MinXY[1] = *min_element(ReProjectionErrorY.begin(), ReProjectionErrorY.end());
		res.MaxXY[0] = *max_element(ReProjectionErrorX.begin(), ReProjectionErrorX.end());
		res.MaxXY[1] = *max_element(ReProjectionErrorY.begin(), ReProjectionErrorY.end());

		res.MeanABSXY[0] = MeanArray(ReProjectionErrorX);
		res.MeanABSXY[1] = MeanArray(ReProjectionErrorY);

		res.STDXY[0] = sqrt(VarianceArray(ReProjectionErrorX, res.MeanABSXY[0]));
		res.STDXY[1] = sqrt(VarianceArray(ReProjectionErrorY, res.MeanABSXY[1]));

		return;
	}
	void VisSfMLens2OPENCV(vector<CameraData> &camera, vector< vector<double> > &xyz)
	{
		const int nCameras = camera.size();

#pragma omp parallel for
		for (int camID = 0; camID < nCameras; camID++)
		{
			ceres::Problem problem;
			CameraData *cam = &camera[camID];
			//if (!cam->available || !cam->convert_1st_radial)
			if (!cam->available)
				continue;

#pragma omp critical
			cout << "\n\t" << "cam: " << cam->filename;


			double
				*focal = cam->FocalLength,
				*imcenter = cam->OpticalCenter,
				*skew = &cam->Skew,
				*rad_1st = &cam->Radialfirst,
				*rad_other = cam->Radialothers,
				*tangential = cam->Tangential,
				*prism = cam->Prism,
				*angleaxis = cam->AngleAxis,
				*trans = cam->Translation;

			const int nPoints = cam->point2D.size();

			for (int i = 0; i < nPoints; i++)
			{
				double observed_x = cam->point2D[i][0];
				double observed_y = cam->point2D[i][1];

				int     ptID = cam->ptID[i];
				double *pt3D = &xyz[ptID][0];

				ceres::CostFunction* cost_function = PinholeReprojectionError::Create(observed_x, observed_y);
				problem.AddResidualBlock(cost_function, NULL, focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D);


				problem.SetParameterBlockConstant(pt3D);
			}
			problem.SetParameterBlockConstant(focal);
			problem.SetParameterBlockConstant(imcenter);
			problem.SetParameterBlockConstant(skew);
			problem.SetParameterBlockConstant(rad_other);
			problem.SetParameterBlockConstant(tangential);
			problem.SetParameterBlockConstant(prism);
			problem.SetParameterBlockConstant(angleaxis);
			problem.SetParameterBlockConstant(trans);


			ceres::Solver::Summary summary;
			ceres::Solver::Options options;

			options.linear_solver_type = ceres::DENSE_SCHUR;
			//options.trust_region_strategy_type = ceres::DOGLEG;
			options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
			options.use_nonmonotonic_steps = true;
			options.minimizer_progress_to_stdout = false;

			ceres::Solve(options, &problem, &summary);

		}
		cout << endl;
	}
	void setCeresOption(const NVM &nvmdata, ceres::Solver::Options &options)
	{
		SYSTEM_INFO sysinfo;
		GetSystemInfo(&sysinfo);
		int nCPUs = sysinfo.dwNumberOfProcessors;

		options.num_threads = nCPUs;
		options.num_linear_solver_threads = nCPUs;
		options.max_num_iterations = 200;
		options.minimizer_progress_to_stdout = true;
		//options.function_tolerance         = 1e-10;
		//options.gradient_tolerance         = 1e-10;

		if (nvmdata.nCamera < 200)
		{
			options.linear_solver_type = ceres::DENSE_SCHUR;
			options.trust_region_strategy_type = ceres::DOGLEG;
			options.use_nonmonotonic_steps = false;
		}
		else
		{
			options.linear_solver_type = ceres::ITERATIVE_SCHUR;
			options.preconditioner_type = ceres::PreconditionerType::CLUSTER_JACOBI;
			//options.visibility_clustering_type = ceres::VisibilityClusteringType::SINGLE_LINKAGE;
		}
	}
	void runBundleAdjustment(vector<CameraData> &camera, vector< vector<double> > &xyz, const ceres::Solver::Options &options, ceres::Solver::Summary &summary, const double thresh)
	{
		const int nCameras = camera.size();

		// problem setting
		ceres::Problem problem;
		for (int camID = 0; camID < nCameras; camID++)
		{
			CameraData *cam = &camera[camID];
			if (!cam->available)
				continue;

			double
				*focal = cam->FocalLength,
				*imcenter = cam->OpticalCenter,
				*skew = &cam->Skew,
				*rad_1st = &cam->Radialfirst,
				*rad_other = cam->Radialothers,
				*tangential = cam->Tangential,
				*prism = cam->Prism,
				*angleaxis = cam->AngleAxis,
				*trans = cam->Translation;

			const int nPoints = cam->point2D.size();

			for (int i = 0; i < nPoints; i++)
			{
				double observed_x = cam->point2D[i][0];
				double observed_y = cam->point2D[i][1];

				int     ptID = cam->ptID[i];
				double *pt3D = &xyz[ptID][0];

				if (thresh > 0)
				{
					double reprojected_pt[2];
					PinholeReprojection(focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D, reprojected_pt);

					double residual_x = reprojected_pt[0] - observed_x;
					double residual_y = reprojected_pt[1] - observed_y;
					double distance = sqrt(residual_x*residual_x + residual_y*residual_y);

					if (distance > thresh)
					{
						cam->inlier[i] = false;
						continue;
					}
				}

				ceres::CostFunction* cost_function = PinholeReprojectionError::Create(observed_x, observed_y);
				problem.AddResidualBlock(cost_function, NULL, focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D);

				setConstantParams(*cam, problem);
			}
		}

		ceres::Solve(options, &problem, &summary);
	}
	void runBundleAdjustmentSharedIntrinsic(vector<CameraData> &camera, vector< vector<double> > &xyz, const ceres::Solver::Options &options, ceres::Solver::Summary &summary, const double thresh)
	{
		const int nCameras = camera.size();

		//Set camera 0 as master camera
		CameraData *cam = &camera[0];
		if (!cam->available)
		{
			printf("Master camera (#0) is not available!\n");
			return;
		}
		double *focal = cam->FocalLength, *imcenter = cam->OpticalCenter, *skew = &cam->Skew,
			*rad_1st = &cam->Radialfirst, *rad_other = cam->Radialothers, *tangential = cam->Tangential, *prism = cam->Prism,
			*omega = &cam->omega, *DistCtr = cam->DistCtr;

		// problem setting
		ceres::Problem problem;
		for (int camID = 0; camID < nCameras; camID++)
		{
			CameraData *cam = &camera[camID];
			if (!cam->available)
				continue;

			double *angleaxis = cam->AngleAxis, *trans = cam->Translation;

			const int nPoints = cam->point2D.size();
			for (int i = 0; i < nPoints; i++)
			{
				double observed_x = cam->point2D[i][0], observed_y = cam->point2D[i][1];

				int     ptID = cam->ptID[i];
				double *pt3D = &xyz[ptID][0];

				if (thresh > 0)
				{
					double reprojected_pt[2];
					if (cam->lenstype == BA_LENS_PINHOLE)
						PinholeReprojection(focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D, reprojected_pt);
					else
						FOVReprojection(focal, imcenter, skew, omega, DistCtr, angleaxis, trans, pt3D, reprojected_pt);

					double residual_x = reprojected_pt[0] - observed_x;
					double residual_y = reprojected_pt[1] - observed_y;
					double distance = sqrt(residual_x*residual_x + residual_y*residual_y);

					if (distance > thresh)
					{
						cam->inlier[i] = false;
						continue;
					}
				}

				if (cam->lenstype == BA_LENS_PINHOLE)
				{
					ceres::CostFunction* cost_function = PinholeReprojectionError::Create(observed_x, observed_y);
					problem.AddResidualBlock(cost_function, NULL, focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D);
				}
				else
				{
					ceres::CostFunction* cost_function = FOVReprojectionError::Create(observed_x, observed_y);
					problem.AddResidualBlock(cost_function, NULL, focal, imcenter, skew, omega, DistCtr, angleaxis, trans, pt3D);
				}
			}
		}

		cam = &camera[0];
		switch (cam->opt_intrinsic)
		{
		case BA_OPT_INTRINSIC_ALL_FIXED:
			problem.SetParameterBlockConstant(focal);
			problem.SetParameterBlockConstant(imcenter);
			problem.SetParameterBlockConstant(skew);
			break;

		case BA_OPT_INTRINSIC_CENTER_FIXED:
			problem.SetParameterBlockConstant(imcenter);
			break;

		case BA_OPT_INTRINSIC_SKEW_ZERO_FIXED:
		case BA_OPT_INTRINSIC_SKEW_FIXED:
			problem.SetParameterBlockConstant(skew);
			break;
		}

		if (cam->lenstype == BA_LENS_PINHOLE)
			switch (cam->opt_lensdistortion)
		{
			case BA_OPT_LENSDIST_RADIAL_AND_TANGENTIAL:
				problem.SetParameterBlockConstant(prism);
				break;

			case BA_OPT_LENSDIST_ALL_FIXED:
				problem.SetParameterBlockConstant(rad_1st);
				problem.SetParameterBlockConstant(rad_other);
				problem.SetParameterBlockConstant(tangential);
				problem.SetParameterBlockConstant(prism);
				break;


			case BA_OPT_LENSDIST_RADIAL_1ST_ONLY:
				problem.SetParameterBlockConstant(rad_other);
				problem.SetParameterBlockConstant(tangential);
				problem.SetParameterBlockConstant(prism);
				break;

			case BA_OPT_LENSDIST_RADIAL_ONLY:
				problem.SetParameterBlockConstant(tangential);
				problem.SetParameterBlockConstant(prism);
				break;
		}

		ceres::Solve(options, &problem, &summary);
	}
	void runBundleAdjustmentPartiallySharedIntrinsic(vector<CameraData> &camera, vector< vector<double> > &xyz, const ceres::Solver::Options &options, ceres::Solver::Summary &summary, const double thresh)
	{
		const int nCameras = camera.size();

		//Set camera 0 as master camera
		CameraData *cam = &camera[0];
		if (!cam->available)
		{
			printf("Master camera (#0) is not available!\n");
			return;
		}
		bool firstSharedCamAdded = false;
		double *focalS, *imcenterS, *skewS,
			*rad_1stS, *rad_otherS, *tangentialS, *prismS,
			*omegaS, *DistCtrS;

		// problem setting
		ceres::Problem problem;
		for (int camID = 0; camID < nCameras; camID++)
		{
			CameraData *cam = &camera[camID];
			if (!cam->available)
				continue;

			if (cam->sharedIntrinisc == true && !firstSharedCamAdded)
				focalS = cam->FocalLength, imcenterS = cam->OpticalCenter, skewS = &cam->Skew,
				rad_1stS = &cam->Radialfirst, rad_otherS = cam->Radialothers, tangentialS = cam->Tangential, prismS = cam->Prism,
				omegaS = &cam->omega, DistCtrS = cam->DistCtr;


			double *focal = cam->FocalLength, *imcenter = cam->OpticalCenter, *skew = &cam->Skew,
				*rad_1st = &cam->Radialfirst, *rad_other = cam->Radialothers, *tangential = cam->Tangential, *prism = cam->Prism,
				*omega = &cam->omega, *DistCtr = cam->DistCtr;

			double *angleaxis = cam->AngleAxis, *trans = cam->Translation;

			const int nPoints = cam->point2D.size();
			for (int i = 0; i < nPoints; i++)
			{
				double observed_x = cam->point2D[i][0], observed_y = cam->point2D[i][1];

				int     ptID = cam->ptID[i];
				double *pt3D = &xyz[ptID][0];

				if (thresh > 0)
				{
					double reprojected_pt[2];
					if (cam->lenstype == BA_LENS_PINHOLE)
						PinholeReprojection(focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D, reprojected_pt);
					else
						FOVReprojection(focal, imcenter, skew, omega, DistCtr, angleaxis, trans, pt3D, reprojected_pt);

					double residual_x = reprojected_pt[0] - observed_x;
					double residual_y = reprojected_pt[1] - observed_y;
					double distance = sqrt(residual_x*residual_x + residual_y*residual_y);

					if (distance > thresh)
					{
						cam->inlier[i] = false;
						continue;
					}
				}

				if (cam->lenstype == BA_LENS_PINHOLE)
				{
					if (cam->sharedIntrinisc == true)
					{
						ceres::CostFunction* cost_function = PinholeReprojectionError::Create(observed_x, observed_y);
						problem.AddResidualBlock(cost_function, NULL, focalS, imcenterS, skewS, rad_1stS, rad_otherS, tangentialS, prismS, angleaxis, trans, pt3D);
					}
					else
					{
						ceres::CostFunction* cost_function = PinholeReprojectionError::Create(observed_x, observed_y);
						problem.AddResidualBlock(cost_function, NULL, focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D);
					}
				}
				else
				{
					if (cam->sharedIntrinisc == true)
					{
						ceres::CostFunction* cost_function = FOVReprojectionError::Create(observed_x, observed_y);
						problem.AddResidualBlock(cost_function, NULL, focalS, imcenterS, skewS, omegaS, DistCtrS, angleaxis, trans, pt3D);
					}
					else
					{
						ceres::CostFunction* cost_function = FOVReprojectionError::Create(observed_x, observed_y);
						problem.AddResidualBlock(cost_function, NULL, focal, imcenter, skew, omega, DistCtr, angleaxis, trans, pt3D);
					}
				}
			}

			//Set up constraints
			if (cam->sharedIntrinisc && !firstSharedCamAdded)
			{
				firstSharedCamAdded = true;
				switch (cam->opt_intrinsic)
				{
				case BA_OPT_INTRINSIC_ALL_FIXED:
					problem.SetParameterBlockConstant(focalS);
					problem.SetParameterBlockConstant(imcenterS);
					problem.SetParameterBlockConstant(skewS);
					break;

				case BA_OPT_INTRINSIC_CENTER_FIXED:
					problem.SetParameterBlockConstant(imcenterS);
					break;

				case BA_OPT_INTRINSIC_SKEW_ZERO_FIXED:
				case BA_OPT_INTRINSIC_SKEW_FIXED:
					problem.SetParameterBlockConstant(skewS);
					break;
				}

				if (cam->lenstype == BA_LENS_PINHOLE)
				{
					switch (cam->opt_lensdistortion)
					{
					case BA_OPT_LENSDIST_RADIAL_AND_TANGENTIAL:
						problem.SetParameterBlockConstant(prismS);
						break;

					case BA_OPT_LENSDIST_ALL_FIXED:
						problem.SetParameterBlockConstant(rad_1stS);
						problem.SetParameterBlockConstant(rad_otherS);
						problem.SetParameterBlockConstant(tangentialS);
						problem.SetParameterBlockConstant(prismS);
						break;


					case BA_OPT_LENSDIST_RADIAL_1ST_ONLY:
						problem.SetParameterBlockConstant(rad_otherS);
						problem.SetParameterBlockConstant(tangentialS);
						problem.SetParameterBlockConstant(prismS);
						break;

					case BA_OPT_LENSDIST_RADIAL_ONLY:
						problem.SetParameterBlockConstant(tangentialS);
						problem.SetParameterBlockConstant(prismS);
						break;
					}
				}
			}
			else if ((cam->sharedIntrinisc && firstSharedCamAdded) && !cam->sharedIntrinisc)
			{
				switch (cam->opt_intrinsic)
				{
				case BA_OPT_INTRINSIC_ALL_FIXED:
					problem.SetParameterBlockConstant(focal);
					problem.SetParameterBlockConstant(imcenter);
					problem.SetParameterBlockConstant(skew);
					break;

				case BA_OPT_INTRINSIC_CENTER_FIXED:
					problem.SetParameterBlockConstant(imcenter);
					break;

				case BA_OPT_INTRINSIC_SKEW_ZERO_FIXED:
				case BA_OPT_INTRINSIC_SKEW_FIXED:
					problem.SetParameterBlockConstant(skew);
					break;
				}

				if (cam->lenstype == BA_LENS_PINHOLE)
				{
					switch (cam->opt_lensdistortion)
					{
					case BA_OPT_LENSDIST_RADIAL_AND_TANGENTIAL:
						problem.SetParameterBlockConstant(prism);
						break;

					case BA_OPT_LENSDIST_ALL_FIXED:
						problem.SetParameterBlockConstant(rad_1st);
						problem.SetParameterBlockConstant(rad_other);
						problem.SetParameterBlockConstant(tangential);
						problem.SetParameterBlockConstant(prism);
						break;


					case BA_OPT_LENSDIST_RADIAL_1ST_ONLY:
						problem.SetParameterBlockConstant(rad_other);
						problem.SetParameterBlockConstant(tangential);
						problem.SetParameterBlockConstant(prism);
						break;

					case BA_OPT_LENSDIST_RADIAL_ONLY:
						problem.SetParameterBlockConstant(tangential);
						problem.SetParameterBlockConstant(prism);
						break;
					}
				}
			}
		}

		ceres::Solve(options, &problem, &summary);
	}
	void runBundleAdjustment(vector<CameraData> &camera, vector< vector<double> > &xyz, vector< vector<bool> > &visMap, const double thresh, const ceres::Solver::Options &options, ceres::Solver::Summary &summary)
	{
		//remove outliers
		if (thresh > 0)
			checkOutlier(camera, xyz, visMap, thresh);

		// problem setting
		const int nCameras = camera.size();
		ceres::Problem problem;
		for (int camID = 0; camID < nCameras; camID++)
		{
			CameraData *cam = &camera[camID];
			if (!cam->available)
				continue;

			double
				*focal = cam->FocalLength,
				*imcenter = cam->OpticalCenter,
				*skew = &cam->Skew,
				*rad_1st = &cam->Radialfirst,
				*rad_other = cam->Radialothers,
				*tangential = cam->Tangential,
				*prism = cam->Prism,
				*angleaxis = cam->AngleAxis,
				*trans = cam->Translation;

			const int nPoints = cam->point2D.size();

			for (int i = 0; i < nPoints; i++)
			{
				int     ptID = cam->ptID[i];
				double *pt3D = &xyz[ptID][0];

				if (!visMap[ptID][camID])
					continue;


				double observed_x = cam->point2D[i][0];
				double observed_y = cam->point2D[i][1];

				ceres::CostFunction* cost_function = PinholeReprojectionError::Create(observed_x, observed_y);
				problem.AddResidualBlock(cost_function, NULL, focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D);

				setConstantParams(*cam, problem);
			}

			//setConstantParams(*cam, problem);


		}
		ceres::Solve(options, &problem, &summary);

	}

	//check outliers and update visMap
	void checkOutlier(vector<CameraData> &camera, vector< vector<double> > &xyz, vector< vector<bool> > &visMap, const double thresh)
	{

		const int nCameras = camera.size();
		for (int camID = 0; camID < nCameras; camID++)
		{
			CameraData *cam = &camera[camID];

			double
				*focal = cam->FocalLength,
				*imcenter = cam->OpticalCenter,
				*skew = &cam->Skew,
				*rad_1st = &cam->Radialfirst,
				*rad_other = cam->Radialothers,
				*tangential = cam->Tangential,
				*prism = cam->Prism,
				*angleaxis = cam->AngleAxis,
				*trans = cam->Translation;

			const int nPoints = cam->point2D.size();

			for (int i = 0; i < nPoints; i++)
			{
				if (!cam->available)
				{
					visMap[i][camID] = false;
					continue;
				}


				int     ptID = cam->ptID[i];
				double *pt3D = &xyz[ptID][0];

				if (!visMap[ptID][camID])
					continue;


				double observed_x = cam->point2D[i][0];
				double observed_y = cam->point2D[i][1];


				double reprojected_pt[2];
				PinholeReprojection(focal, imcenter, skew, rad_1st, rad_other, tangential, prism, angleaxis, trans, pt3D, reprojected_pt);

				double residual_x = reprojected_pt[0] - observed_x;
				double residual_y = reprojected_pt[1] - observed_y;
				double distance = sqrt(residual_x*residual_x + residual_y*residual_y);

				if (distance > thresh)
					visMap[ptID][camID] = false;

			}
		}


		//remove a point which is viewed by only 1 camera
		int nPoints = visMap.size();
		for (int ptID = 0; ptID < nPoints; ptID++)
		{
			int nViews = 0;
			for (int camID = 0; camID < nCameras; camID++)
			{
				if (visMap[ptID][camID])
					nViews++;

				if (nViews > 2)
					break;
			}

			if (nViews < 2)
				std::fill(visMap[ptID].begin(), visMap[ptID].end(), false);
			//for (int camID = 0; camID < nCameras; camID++)
			//	visMap[ptID][camID] = false;
		}


	}
	void setConstantParams(CameraData &cam, ceres::Problem &problem)
	{

		double
			*focal = cam.FocalLength,
			*center = cam.OpticalCenter,
			*skew = &cam.Skew,
			*rad_1st = &cam.Radialfirst,
			*rad_other = cam.Radialothers,
			*tangential = cam.Tangential,
			*prism = cam.Prism,
			*angleaxis = cam.AngleAxis,
			*trans = cam.Translation;


		switch (cam.opt_intrinsic)
		{
		case BA_OPT_INTRINSIC_ALL:
			break;

		case BA_OPT_INTRINSIC_ALL_FIXED:
			problem.SetParameterBlockConstant(focal);
			problem.SetParameterBlockConstant(center);
			problem.SetParameterBlockConstant(skew);
			break;

		case BA_OPT_INTRINSIC_CENTER_FIXED:
			problem.SetParameterBlockConstant(center);
			break;

		case BA_OPT_INTRINSIC_SKEW_ZERO_FIXED:
		case BA_OPT_INTRINSIC_SKEW_FIXED:
			problem.SetParameterBlockConstant(skew);
			break;
		}



		switch (cam.opt_lensdistortion)
		{
		case BA_OPT_LENSDIST_RADIAL_AND_TANGENTIAL:
			problem.SetParameterBlockConstant(prism);
			break;

		case BA_OPT_LENSDIST_ALL_FIXED:
			problem.SetParameterBlockConstant(rad_1st);
			problem.SetParameterBlockConstant(rad_other);
			problem.SetParameterBlockConstant(tangential);
			problem.SetParameterBlockConstant(prism);
			break;


		case BA_OPT_LENSDIST_RADIAL_1ST_ONLY:
			problem.SetParameterBlockConstant(rad_other);
			problem.SetParameterBlockConstant(tangential);
			problem.SetParameterBlockConstant(prism);
			break;

		case BA_OPT_LENSDIST_RADIAL_ONLY:
			problem.SetParameterBlockConstant(tangential);
			problem.SetParameterBlockConstant(prism);
			break;
		}


		switch (cam.opt_extrinsic)
		{
		case BA_OPT_EXTRINSIC_ALL:
			break;

		case BA_OPT_EXTRINSIC_ALL_FIXED:
			problem.SetParameterBlockConstant(angleaxis);
			problem.SetParameterBlockConstant(trans);
			break;

		case BA_OPT_EXTRINSIC_R_FIXED:
			problem.SetParameterBlockConstant(angleaxis);
			break;

		case BA_OPT_EXTRINSIC_T_FIXED:
			problem.SetParameterBlockConstant(trans);
			break;
		}

	}
	bool loadNVM(const string filepath, NVM &nvmdata, int sharedIntrinsics)
	{
		ifstream ifs(filepath);
		if (ifs.fail())
		{
			cerr << "Cannot load " << filepath << endl;
			return false;
		}

		string token;
		ifs >> token; //NVM_V3
		if (token != "NVM_V3")
		{
			cerr << "Can only load NVM_V3" << endl;
			return false;
		}

		double fx, fy, u0, v0, radial1;
		if (sharedIntrinsics == 1)
			ifs >> token >> fx >> u0 >> fy >> v0 >> radial1;

		//loading camera parameters
		int nCameras;
		ifs >> nCameras;
		if (nCameras <= 1)
		{
			cerr << "# of cameras must be more than 1." << endl;
			return false;
		}
		nvmdata.nCamera = nCameras;
		for (int camID = 0; camID < nCameras; camID++)
		{
			string filename;
			double f;
			vector<double> q(4), c(3), d(2);
			ifs >> filename >> f >> q[0] >> q[1] >> q[2] >> q[3] >> c[0] >> c[1] >> c[2] >> d[0] >> d[1];

			if (sharedIntrinsics == 1)
				nvmdata.focallength.push_back(0.5*(fx + fy));
			else
				nvmdata.focallength.push_back(f);

			nvmdata.filenames.push_back(filename);
			nvmdata.quaternion.push_back(q);
			nvmdata.position.push_back(c);
			nvmdata.firstradial.push_back(d[0]);

			nvmdata.filename_id[filename] = camID;
		}


		//loading 2D and 3D points
		int nPoints;
		ifs >> nPoints;
		if (nCameras <= 0)
		{
			cerr << "# of 3D points is 0." << endl;
			return false;
		}

		nvmdata.n3dPoints = nPoints;
		for (int i = 0; i < nPoints; i++)
		{
			vector<double> xyz(3);
			vector<int> rgb(3);
			string info;
			ifs >> xyz[0] >> xyz[1] >> xyz[2]
				>> rgb[0] >> rgb[1] >> rgb[2];

			getline(ifs, info);

			nvmdata.xyz.push_back(xyz);
			nvmdata.rgb.push_back(rgb);
			nvmdata.measurementinfo.push_back(info);
		}

		return true;
	}
	bool checkOptionIntrinsics(const int option)
	{
		bool ok = true;

		switch (option)
		{
		case BA_OPT_INTRINSIC_ALL:
		case BA_OPT_INTRINSIC_ALL_FIXED:
		case BA_OPT_INTRINSIC_SKEW_FIXED:
		case BA_OPT_INTRINSIC_SKEW_ZERO_FIXED:
		case BA_OPT_INTRINSIC_CENTER_FIXED:
		case BA_OPT_LOAD_INTRINSICS_BUT_JUST_AS_INITIAL_GUESS:
			break;

		default:
			ok = false;
			break;
		}

		return ok;
	}
	bool checkOptionExtrinsics(const int option)
	{
		bool ok = true;

		switch (option)
		{
		case BA_OPT_EXTRINSIC_ALL:
		case BA_OPT_EXTRINSIC_ALL_FIXED:
		case BA_OPT_EXTRINSIC_R_FIXED:
		case BA_OPT_EXTRINSIC_T_FIXED:
			break;

		default:
			ok = false;
			break;
		}

		return ok;
	}
	bool checkOptionLensDist(const int option)
	{
		bool ok = true;

		switch (option)
		{
		case BA_LENS_FISHEYE:
		case BA_OPT_LENSDIST_RADIAL_AND_TANGENTIAL:
		case BA_OPT_LENSDIST_ALL_FIXED:
		case BA_OPT_LENSDIST_RADIAL_1ST_ONLY:
		case BA_OPT_LENSDIST_RADIAL_ONLY:
			break;

		default:
			ok = false;
			break;
		}

		return ok;
	}
	bool initCameraData(const NVM &nvmdata, const string filepath, vector<CameraData> &camera, int sharedIntrinisc)
	{
		//load image info
		ifstream ifs(filepath);
		if (ifs.fail())
		{
			cerr << "Cannot open " << filepath << endl;
			return false;
		}

		vector<int> imgWidth, imgHeight, opt_intrinsic, opt_lensdistortion, opt_extrinsic;
		vector<bool> available;
		vector<bool> intrinsicShared;
		vector<string> filenames;
		string str;
		char tfile[200];
		if (sharedIntrinisc == 1)
		{
			int w, h, option1, option2, option3;
			getline(ifs, str);
			stringstream ss(str);
			ss >> w >> h >> option1 >> option2 >> option3;

			for (int ii = 0; ii < nvmdata.nCamera; ii++)
			{
				sprintf(tfile, "%d.ppm", ii);
				std::string ifile = tfile;
				if (!checkOptionIntrinsics(option1))
				{
					cerr << "\n" << "Option for intrinsic parameters is not valid in " << ifile << endl;
					return false;
				}
				if (!checkOptionLensDist(option2))
				{
					cerr << "\n" << "Option for lens distortion is not valid in " << ifile << endl;
					return false;
				}
				if (!checkOptionExtrinsics(option3))
				{
					cerr << "\n" << "Option for extrinsic parameters is not valid in " << ifile << endl;
					return false;
				}

				filenames.push_back(ifile);
				imgWidth.push_back(w);
				imgHeight.push_back(h);
				opt_intrinsic.push_back(option1);
				opt_lensdistortion.push_back(option2);
				opt_extrinsic.push_back(option3);

				available.push_back(true);
			}
		}
		else
		{
			while (getline(ifs, str))
			{
				string file;
				int w, h, option1, option2, option3, avail, shared;

				stringstream ss(str);
				ss >> file >> w >> h >> option1 >> option2 >> option3 >> avail >> shared;

				if (!checkOptionIntrinsics(option1))
				{
					cerr << "\n" << "Option for intrinsic parameters is not valid in " << file << endl;
					return false;
				}
				if (!checkOptionLensDist(option2))
				{
					cerr << "\n" << "Option for lens distortion is not valid in " << file << endl;
					return false;
				}
				if (!checkOptionExtrinsics(option3))
				{
					cerr << "\n" << "Option for extrinsic parameters is not valid in " << file << endl;
					return false;
				}

				filenames.push_back(file);
				imgWidth.push_back(w);
				imgHeight.push_back(h);
				opt_intrinsic.push_back(option1);
				opt_lensdistortion.push_back(option2);
				opt_extrinsic.push_back(option3);

				bool flag = avail > 0 ? true : false;
				available.push_back(flag);

				if (shared == 1)
					intrinsicShared.push_back(true);
				else
					intrinsicShared.push_back(false);
			}
		}

		int nCameras = nvmdata.nCamera;
		if (filenames.size() > nCameras)
			cerr << "\n" << "Warning: " << "The number of cameras in the ini file is more than that in the NVM file." << "The ini file may not be valid." << endl;
		else if (filenames.size() < nCameras)
		{
			cerr << "\n" << "Error: " << "The number of cameras in the ini file is less than that in the NVM file.\n" << "The ini file is not valid." << endl;
			return false;
		}

		//initialize cmaera data
		camera.resize(nCameras);
		int nImages = filenames.size();
		int found = 0;
		vector<bool> load_ok(nCameras, false);

		bool first_warning = true;
		for (int i = 0; i < nImages; i++)
		{
			map<string, int>::const_iterator itr = nvmdata.filename_id.find(filenames[i]);
			if (itr == nvmdata.filename_id.end())
			{
				if (first_warning)
				{
					cerr << "\n";
					first_warning = false;
				}

				cerr << "Warning: " << "Cannot find an image file (" << filenames[i] << ") in the NVM file." << endl;
				continue;
			}

			//string imgfile = itr->first; <- is equal to filenames[i]
			int camID = itr->second;

			CameraData *cam = &camera[camID];
			found++;
			load_ok[camID] = true;

			cam->filename = filenames[i];
			cam->imgWidth = imgWidth[i];
			cam->imgHeight = imgHeight[i];

			cam->opt_intrinsic = opt_intrinsic[i];
			cam->opt_lensdistortion = opt_lensdistortion[i];
			cam->opt_extrinsic = opt_extrinsic[i];
			cam->available = available[i];
			cam->sharedIntrinisc = intrinsicShared[i];

			cam->FocalLength[0] = nvmdata.focallength[camID];
			cam->FocalLength[1] = nvmdata.focallength[camID];
			cam->OpticalCenter[0] = 0.5*(cam->imgWidth - 1.0);
			cam->OpticalCenter[1] = 0.5*(cam->imgHeight - 1.0);
			cam->Skew = 0.0;

			if (cam->opt_lensdistortion == -1)
			{
				cam->lenstype = BA_LENS_FISHEYE;
				cam->omega = 0.001, cam->DistCtr[0] = cam->OpticalCenter[0], cam->DistCtr[1] = cam->OpticalCenter[1];
			}
			else
			{
				cam->lenstype = BA_LENS_PINHOLE;
				cam->Radialfirst = -1.0*nvmdata.firstradial[camID]; //good approximation
				for (int j = 0; j < 2; j++)
					cam->Radialothers[j] = 0.0, cam->Tangential[j] = 0.0, cam->Prism[j] = 0.0;
			}

			//quaternion to angle-axis
			double q[4] = { nvmdata.quaternion[camID][0], nvmdata.quaternion[camID][1], nvmdata.quaternion[camID][2], nvmdata.quaternion[camID][3] };
			double angle[3];
			ceres::QuaternionToAngleAxis(q, angle);
			cam->AngleAxis[0] = angle[0];
			cam->AngleAxis[1] = angle[1];
			cam->AngleAxis[2] = angle[2];

			//position to translation t=-R*c
			double c[3] = { nvmdata.position[camID][0], nvmdata.position[camID][1], nvmdata.position[camID][2] };
			double t[3];
			ceres::QuaternionRotatePoint(q, c, t);

			cam->Translation[0] = -t[0];
			cam->Translation[1] = -t[1];
			cam->Translation[2] = -t[2];
		}

		if (found != nCameras)
		{
			cerr << "Error: " << "The ini file does not have information for following cameras:" << endl;
			for (int camID = 0; camID < nCameras; camID++)
			{
				if (!load_ok[camID])
					cerr << nvmdata.filenames[camID] << endl;
			}

			return false;
		}

		//copy point data to CameraData
		for (int ptID = 0; ptID < nvmdata.n3dPoints; ptID++)
		{
			stringstream info(nvmdata.measurementinfo[ptID]);
			int nObserbed;

			info >> nObserbed;
			for (int i = 0; i < nObserbed; i++)
			{
				int camID, fID;
				double observed_x, observed_y;
				info >> camID >> fID >> observed_x >> observed_y;

				CameraData *cam = &camera[camID];
				if (!cam->available)
					continue;


				vector<double> point2D(2);
				point2D[0] = observed_x + 0.5*(cam->imgWidth - 1.0);
				point2D[1] = observed_y + 0.5*(cam->imgHeight - 1.0);
				cam->point2D.push_back(point2D);

				cam->ptID.push_back(ptID);
				//cam->fID.push_back(fID);

				cam->inlier.push_back(true);
			}

		}

		return true;
	}
	bool loadInitialIntrinsics(const string intrinsicfile, const map<string, int> &filename_id, vector<CameraData> &camera, int sharedIntrinsics)
	{
		ifstream ifs(intrinsicfile);
		if (ifs.fail())
		{
			cerr << "Cannot open " << intrinsicfile << endl;
			return false;
		}


		if (sharedIntrinsics == 1)
		{
			bool first_warning = true;
			string str;
			while (getline(ifs, str))
			{
				double fx, fy, s, u0, v0, a0, a1, a2, p0, p1, s0, s1;
				stringstream ss(str);
				ss >> fx >> fy >> s >> u0 >> v0
					>> a0 >> a1 >> a2
					>> p0 >> p1 >> s0 >> s1;

				for (int camID = 0; camID < camera.size(); camID++)
				{
					camera[camID].FocalLength[0] = fx, camera[camID].FocalLength[1] = fy, camera[camID].Skew = s, camera[camID].OpticalCenter[0] = u0, camera[camID].OpticalCenter[1] = v0;

					camera[camID].Radialfirst = a0, camera[camID].Radialothers[0] = a1, camera[camID].Radialothers[1] = a2;
					camera[camID].Tangential[0] = p0, camera[camID].Tangential[1] = p0;
					camera[camID].Prism[0] = s0, camera[camID].Prism[1] = s1;


					if (camera[camID].opt_intrinsic == BA_OPT_LOAD_INTRINSICS_BUT_JUST_AS_INITIAL_GUESS)
						camera[camID].opt_intrinsic = BA_OPT_INTRINSIC_ALL;
					else
					{
						camera[camID].opt_intrinsic = BA_OPT_INTRINSIC_ALL_FIXED;
						camera[camID].opt_lensdistortion = BA_OPT_LENSDIST_ALL_FIXED;
					}
				}
			}
		}
		else
		{
			bool first_warning = true;
			string str;
			while (getline(ifs, str))
			{
				string file;
				int lenstype, w, h;
				double fx, fy, s, u0, v0, a0, a1, a2, p0, p1, s0, s1;
				stringstream ss(str);
				ss >> file >> lenstype>>w>>h
					>> fx >> fy >> s >> u0 >> v0
					>> a0 >> a1 >> a2
					>> p0 >> p1 >> s0 >> s1;

				map<string, int>::const_iterator itr = filename_id.find(file);
				if (itr == filename_id.end())
				{
					if (first_warning)
					{
						cerr << "\n";
						first_warning = false;
					}

					cerr << "Warning: " << "Cannot find an image file (" << file << ") in the NVM file. Ignore it and continue." << endl;
					continue;
				}

				int camID = itr->second;

				camera[camID].lenstype = lenstype;
				camera[camID].FocalLength[0] = fx, camera[camID].FocalLength[1] = fy, camera[camID].Skew = s, camera[camID].OpticalCenter[0] = u0, camera[camID].OpticalCenter[1] = v0;

				camera[camID].Radialfirst = a0, camera[camID].Radialothers[0] = a1, camera[camID].Radialothers[1] = a2;
				camera[camID].Tangential[0] = p0, camera[camID].Tangential[1] = p0;
				camera[camID].Prism[0] = s0, camera[camID].Prism[1] = s1;


				if (camera[camID].opt_intrinsic == BA_OPT_LOAD_INTRINSICS_BUT_JUST_AS_INITIAL_GUESS)
					camera[camID].opt_intrinsic = BA_OPT_INTRINSIC_ALL;
				else
				{
					camera[camID].opt_intrinsic = BA_OPT_INTRINSIC_ALL_FIXED;
					camera[camID].opt_lensdistortion = BA_OPT_LENSDIST_ALL_FIXED;
				}
			}
		}

		return true;
	}
	bool loadIntrinsics(const string intrinsicfile, vector<CameraData> &camera)
	{
		ifstream ifs(intrinsicfile);
		if (ifs.fail())
		{
			cerr << "Cannot open " << intrinsicfile << endl;
			return false;
		}

		bool empty = camera.empty() ? true : false;

		int camID = 0;
		string str;
		while (getline(ifs, str))
		{
			string filename;
			double fx, fy, s, u0, v0, a0, a1, a2, p0, p1, s0, s1;
			stringstream ss(str);
			ss >> filename
				>> fx >> fy >> s >> u0 >> v0
				>> a0 >> a1 >> a2
				>> p0 >> p1 >> s0 >> s1;

			CameraData cam;

			cam.filename = filename;
			cam.FocalLength[0] = fx;
			cam.FocalLength[1] = fy;
			cam.OpticalCenter[0] = u0;
			cam.OpticalCenter[1] = v0;
			cam.Skew = s;

			cam.Radialfirst = a0;
			cam.Radialothers[0] = a1;
			cam.Radialothers[1] = a2;
			cam.Tangential[0] = p0;
			cam.Tangential[1] = p0;
			cam.Prism[0] = s0;
			cam.Prism[1] = s1;

			cam.available = true;

			if (empty)
				camera.push_back(cam);
			else
				camera.at(camID) = cam;


			camID++;
		}


		return true;
	}
	bool loadExtrinsics(const string extrinsicfile, vector<CameraData> &camera)
	{
		ifstream ifs(extrinsicfile);
		if (ifs.fail())
		{
			cerr << "Cannot open " << extrinsicfile << endl;
			return false;
		}

		bool empty = camera.empty() ? true : false;

		int camID = 0;
		string str;
		while (getline(ifs, str))
		{
			string filename;
			double r0, r1, r2, t0, t1, t2;
			stringstream ss(str);
			ss >> filename
				>> r0 >> r1 >> r2
				>> t0 >> t1 >> t2;


			CameraData cam;

			cam.AngleAxis[0] = r0;
			cam.AngleAxis[1] = r1;
			cam.AngleAxis[2] = r2;
			cam.Translation[0] = t0;
			cam.Translation[1] = t1;
			cam.Translation[2] = t2;

			cam.available = true;

			if (empty)
				camera.push_back(cam);
			else
				camera.at(camID) = cam;


			camID++;
		}


		return true;
	}
	bool loadAllCameraParams(const string cameraparamfile, vector<CameraData> &camera)
	{
		ifstream ifs(cameraparamfile);
		if (ifs.fail())
		{
			cerr << "Cannot open " << cameraparamfile << endl;
			return false;
		}


		bool empty = camera.empty() ? true : false;

		int camID = 0;
		string str;
		while (getline(ifs, str))
		{
			string cameraname;
			double
				fx, fy, s, u0, v0, // focal length (x, y), skew, optical center (x, y)
				a0, a1, a2,        // radial distortion
				p0, p1,            // tangential distortion
				s0, s1,            // prism distortion
				r0, r1, r2,        // angle-axis rotation parameters
				t0, t1, t2;        // translation


			stringstream ss(str);
			ss >> cameraname
				>> fx >> fy >> s >> u0 >> v0
				>> a0 >> a1 >> a2
				>> p0 >> p1
				>> s0 >> s1
				>> r0 >> r1 >> r2
				>> t0 >> t1 >> t2;

			BA::CameraData cam;
			cam.filename = cameraname;
			cam.FocalLength[0] = fx;
			cam.FocalLength[1] = fy;
			cam.Skew = s;
			cam.OpticalCenter[0] = u0;
			cam.OpticalCenter[1] = v0;

			cam.Radialfirst = a0;
			cam.Radialothers[0] = a1;
			cam.Radialothers[1] = a2;
			cam.Tangential[0] = p0;
			cam.Tangential[1] = p1;
			cam.Prism[0] = s0;
			cam.Prism[1] = s1;

			cam.AngleAxis[0] = r0;
			cam.AngleAxis[1] = r1;
			cam.AngleAxis[2] = r2;
			cam.Translation[0] = t0;
			cam.Translation[1] = t1;
			cam.Translation[2] = t2;

			if (empty)
				camera.push_back(cam);
			else
				camera.at(camID) = cam;

			camID++;
		}

		return true;

	}
	bool saveCameraAllParams(const string filename, const string sep, const vector<CameraData> &camera)
	{
		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}

		ofs << scientific << setprecision(16);

		int notValidCount = 0;
		for (int ii = 0; ii < camera.size(); ii++)
			if (!camera[ii].available)
				notValidCount++;

		ofs << camera.size() - notValidCount << endl;
		for (int camID = 0; camID < camera.size(); camID++)
		{
			if (!camera[camID].available)
				continue;

string file = camera[camID].filename;

double fx = camera[camID].FocalLength[0];
double fy = camera[camID].FocalLength[1];
double s = camera[camID].Skew;
double u0 = camera[camID].OpticalCenter[0];
double v0 = camera[camID].OpticalCenter[1];

double r0 = camera[camID].AngleAxis[0];
double r1 = camera[camID].AngleAxis[1];
double r2 = camera[camID].AngleAxis[2];
double t0 = camera[camID].Translation[0];
double t1 = camera[camID].Translation[1];
double t2 = camera[camID].Translation[2];

if (camera[camID].lenstype == BA_LENS_PINHOLE)
{
	double a0 = camera[camID].Radialfirst;
	double a1 = camera[camID].Radialothers[0];
	double a2 = camera[camID].Radialothers[1];
	double p0 = camera[camID].Tangential[0];
	double p1 = camera[camID].Tangential[1];
	double s0 = camera[camID].Prism[0];
	double s1 = camera[camID].Prism[1];
	ofs << file << sep << camera[0].lenstype << sep << camera[camID].imgWidth << sep << camera[camID].imgHeight << sep
		<< fx << sep << fy << sep << s << sep << u0 << sep << v0 << sep
		<< a0 << sep << a1 << sep << a2 << sep
		<< p0 << sep << p1 << sep << s0 << sep << s1 << sep
		<< r0 << sep << r1 << sep << r2 << sep
		<< t0 << sep << t1 << sep << t2 << endl;
}
else
{
	double omega = camera[camID].omega, DistCtrX = camera[camID].DistCtr[0], DistCtrY = camera[camID].DistCtr[1];
	ofs << file << sep << camera[0].lenstype << sep << camera[camID].imgWidth << sep << camera[camID].imgHeight << sep
		<< fx << sep << fy << sep << s << sep << u0 << sep << v0 << sep
		<< omega << sep << DistCtrX << sep << DistCtrY << sep
		<< r0 << sep << r1 << sep << r2 << sep
		<< t0 << sep << t1 << sep << t2 << endl;
}
		}

		return true;
	}
	bool saveCameraExtrinsics(const string filename, const string sep, const vector<CameraData> &camera)
	{
		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}

		ofs << scientific << setprecision(16);
		for (int camID = 0; camID < camera.size(); camID++)
		{
			if (!camera[camID].available)
				continue;

			string file = camera[camID].filename;

			double r0 = camera[camID].AngleAxis[0];
			double r1 = camera[camID].AngleAxis[1];
			double r2 = camera[camID].AngleAxis[2];
			double t0 = camera[camID].Translation[0];
			double t1 = camera[camID].Translation[1];
			double t2 = camera[camID].Translation[2];

			ofs << file << sep
				<< r0 << sep << r1 << sep << r2 << sep
				<< t0 << sep << t1 << sep << t2 << endl;
		}

		return true;
	}

	bool SaveIntrinsicResults(const string filename, const string sep, const vector<CameraData> &AllViewsParas)
	{
		//Note that visCamualSfm use different lens model than openCV or matlab or yours (inverse model)
		int id = 0, LensType;
		double fx, fy, skew, u0, v0, r0, r1, r2, t0, t1, omega, DistCtrX, DistCtrY;

		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}
		ofs << scientific << setprecision(16);

		int nCams = AllViewsParas.size();
		for (int ii = 0; ii < nCams; ii++)
		{
			int id = ii;
			for (int jj = 0; nCams; jj++)
			{
				string filename = AllViewsParas[jj].filename;
				std::size_t pos = AllViewsParas[jj].filename.find(".ppm");
				filename.erase(pos, 4);
				const char * str = filename.c_str();
				int viewId = atoi(str);
				if (viewId == id)
				{
					id = jj;
					break;
				}
			}

			LensType = AllViewsParas[id].lenstype;
			fx = AllViewsParas[id].FocalLength[0], fy = AllViewsParas[id].FocalLength[1], skew = AllViewsParas[id].Skew, u0 = AllViewsParas[id].OpticalCenter[0], v0 = AllViewsParas[id].OpticalCenter[1];

			if (LensType == BA_LENS_PINHOLE)
			{
				r0 = AllViewsParas[id].Radialfirst, r1 = AllViewsParas[id].Radialothers[0], r2 = AllViewsParas[id].Radialothers[1];
				t0 = AllViewsParas[id].Tangential[0], t1 = AllViewsParas[id].Tangential[1];

				ofs << LensType << sep << AllViewsParas[id].imgWidth << sep << AllViewsParas[id].imgHeight << sep
					<< fx << sep << fy << sep << skew << sep << u0 << sep << v0 << sep
					<< r0 << sep << r1 << sep << r2 << sep
					<< t0 << sep << t1 << sep << 0.0 << sep << 0.0 << endl;
			}
			else
			{
				omega = AllViewsParas[id].omega, DistCtrX = AllViewsParas[id].DistCtr[0], DistCtrY = AllViewsParas[id].DistCtr[1];

				ofs << LensType << sep << AllViewsParas[id].imgWidth << sep << AllViewsParas[id].imgHeight << sep
					<< fx << sep << fy << sep << skew << sep << u0 << sep << v0 << sep
					<< omega << sep << DistCtrX << sep << DistCtrY << endl;
			}
		}

		return true;
	}
	bool saveCameraIntrinsics(const string filename, const string sep, const vector<CameraData> &camera)
	{
		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}

		ofs << scientific << setprecision(16);
		for (int camID = 0; camID < camera.size(); camID++)
		{
			if (!camera[camID].available)
				continue;

			string file = camera[camID].filename;

			double fx = camera[camID].FocalLength[0];
			double fy = camera[camID].FocalLength[1];
			double s = camera[camID].Skew;
			double u0 = camera[camID].OpticalCenter[0];
			double v0 = camera[camID].OpticalCenter[1];

			if (camera[camID].lenstype == BA_LENS_PINHOLE)
			{
				double a0 = camera[camID].Radialfirst;
				double a1 = camera[camID].Radialothers[0];
				double a2 = camera[camID].Radialothers[1];
				double p0 = camera[camID].Tangential[0];
				double p1 = camera[camID].Tangential[1];
				double s0 = camera[camID].Prism[0];
				double s1 = camera[camID].Prism[1];
				ofs << file << sep << camera[camID].lenstype << sep
					<< fx << sep << fy << sep << s << sep << u0 << sep << v0 << sep
					<< a0 << sep << a1 << sep << a2 << sep
					<< p0 << sep << p1 << sep << s0 << sep << s1 << endl;
			}
			else
			{
				double omega = camera[camID].omega, DistCtrX = camera[camID].DistCtr[0], DistCtrY = camera[camID].DistCtr[1];
				ofs << file << sep << camera[camID].lenstype << sep
					<< fx << sep << fy << sep << s << sep << u0 << sep << v0 << sep
					<< omega << sep << DistCtrX << sep << DistCtrY << endl;
			}
		}

		return true;
	}
	bool save3Dpoints(const string filename, const string sep, const vector< vector<double> > &xyz)
	{
		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}
		ofs << scientific << setprecision(16);
		for (int ptID = 0; ptID < xyz.size(); ptID++)
		{
			double x = xyz[ptID][0];
			double y = xyz[ptID][1];
			double z = xyz[ptID][2];

			ofs << x << sep << y << sep << z << endl;
		}

		return true;
	}
	bool saveReprojectionError(const string filename, const string sep, const residualData &res, const vector<CameraData> &camera, const int order)
	{
		ofstream ofs(filename);
		ofs << scientific << setprecision(5);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}
		for (int i = 0; i < res.error.size(); i++)
		{
			int    ptID = res.ID[i][0];
			int    camID = res.ID[i][1];

			string cameraname = camera[camID].filename;

			double x0 = res.observed_pt[i][0];
			double y0 = res.observed_pt[i][1];

			double x1 = res.reprojected_pt[i][0];
			double y1 = res.reprojected_pt[i][1];

			double res_x = res.error[i][0];
			double res_y = res.error[i][1];

			if (order == 0)
				ofs << cameraname << sep << ptID << sep;
			else
				ofs << ptID << sep << cameraname << sep;


			ofs << x0 << sep << y0 << sep
				<< x1 << sep << y1 << sep
				<< res_x << sep << res_y << endl;
		}

		return true;
	}
	bool saveAllData(const string path, const vector<CameraData> &camera, const vector< vector<double> > &xyz, const residualData &res, const string prefix, const bool After)
	{
		string postfix;
		if (After)
			postfix = "_after";
		else
			postfix = "_before";


		string sep = "\t";  //separator
		string ext = ".txt";//extension

		string filename;


		filename = path + "/" + prefix + "Camera_AllParams" + postfix + ext;
		saveCameraAllParams(filename, sep, camera);


		//filename = path + "/" + prefix + "Camera_Extrinsics" + postfix + ext;
		//saveCameraExtrinsics(filename, sep, camera);


		filename = path + "/" + prefix + "Camera_Intrinsics" + postfix + ext;
		saveCameraIntrinsics(filename, sep, camera);

		//filename = path + "/" + "DevicesIntrinsics" +  ext;
		//SaveIntrinsicResults(filename, sep, camera);

		filename = path + "/" + prefix + "ReprojectionError" + postfix + ext;
		saveReprojectionError(filename, sep, res, camera, 0);


		return true;
	}
	void dispCameraParams(const CameraData &camera)
	{
		cout << "\n"
			<< "Name: " << camera.filename << "\n"
			<< "Available   : " << camera.available << "\n"
			<< "Focal Length: " << camera.FocalLength[0] << ", " << camera.FocalLength[1] << "\n"
			<< "Image Center: " << camera.OpticalCenter[0] << ", " << camera.OpticalCenter[1] << "\n"
			<< "Skew        : " << camera.Skew << "\n"
			<< "Radial      : " << camera.Radialfirst << ", " << camera.Radialothers[0] << ", " << camera.Radialothers[1] << ", " << "\n"
			<< "Tangential  : " << camera.Tangential[0] << ", " << camera.Tangential[1] << "\n"
			<< "Prism       : " << camera.Prism[0] << ", " << camera.Prism[1] << "\n"
			<< "\n"
			<< "Angle Axis  : " << camera.AngleAxis[0] << ", " << camera.AngleAxis[1] << ", " << camera.AngleAxis[2] << "\n"
			<< "Translation : " << camera.Translation[0] << ", " << camera.Translation[1] << ", " << camera.Translation[2] << "\n" << endl;
	}
	bool saveNVM(const string path, const string inputnvmname, const vector<CameraData> &camera, const vector< vector<double> > &xyz, const NVM &nvmdata)
	{
		string filename = path + "/BA_" + inputnvmname;
		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}


		// camera information
		//int nCameras = camera.size();
		int nCameras = 0;
		for (int camID = 0; camID < camera.size(); camID++)
			if (camera[camID].available)
				nCameras++;


		ofs << "NVM_V3 \n"
			<< "\n"
			<< nCameras << endl;

		ofs << fixed << setprecision(12);

		for (int camID = 0; camID < nCameras; camID++)
		{
			const CameraData *cam = &camera[camID];
			if (!cam->available)
				continue;

			string file = cam->filename;
			double fx = cam->FocalLength[0];
			double rad = -cam->Radialfirst;


			// angle axis -> quaternion
			double angle[3] = { cam->AngleAxis[0], cam->AngleAxis[1], cam->AngleAxis[2] };
			double q[4];
			ceres::AngleAxisToQuaternion(angle, q);


			// translation -> camera position c=-R'*t
			double t[3] = { cam->Translation[0], cam->Translation[1], cam->Translation[2] };
			double c[3];
			angle[0] = -angle[0];
			angle[1] = -angle[1];
			angle[2] = -angle[2];
			ceres::AngleAxisRotatePoint(angle, t, c);
			c[0] = -c[0];
			c[1] = -c[1];
			c[2] = -c[2];


			ofs << file << "\t"
				<< fx << " "
				<< q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " "
				<< c[0] << " " << c[1] << " " << c[2] << " "
				<< rad << " " << 0 << "\n";
		}
		ofs << endl;


		// point information
		const int nPoints = xyz.size();
		ofs << nPoints << endl;

		for (int ptID = 0; ptID < nPoints; ptID++)
		{
			double x = xyz[ptID][0];
			double y = xyz[ptID][1];
			double z = xyz[ptID][2];

			int r = nvmdata.rgb[ptID][0];
			int g = nvmdata.rgb[ptID][1];
			int b = nvmdata.rgb[ptID][2];

			string info = nvmdata.measurementinfo[ptID];

			ofs << x << " " << y << " " << z << " "
				<< r << " " << g << " " << b
				<< info
				<< endl;
		}
		ofs << "\n"
			<< "\n"
			<< "\n"
			<< "0\n"
			<< "\n"
			<< "#the last part of NVM file points to the PLY files\n"
			<< "#the first number is the number of associated PLY files\n"
			<< "#each following number gives a model - index that has PLY\n"
			<< "0"
			<< endl;

		return true;
	}
	bool saveNVM(const string path, const string outputnvmname, const vector<CameraData> &camera, const vector< vector<double> > &xyz)
	{
		string filename = path + "/" + outputnvmname;
		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}


		// camera information
		//int nCameras = camera.size();
		int nCameras = 0;
		for (int camID = 0; camID < camera.size(); camID++)
			if (camera[camID].available)
				nCameras++;


		ofs << "NVM_V3 \n"
			<< "\n"
			<< nCameras << endl;

		ofs << fixed << setprecision(12);

		for (int camID = 0; camID < nCameras; camID++)
		{
			const CameraData *cam = &camera[camID];
			if (!cam->available)
				continue;

			string file = cam->filename;
			double fx = cam->FocalLength[0];
			double rad = -cam->Radialfirst;


			// angle axis -> quaternion
			double angle[3] = { cam->AngleAxis[0], cam->AngleAxis[1], cam->AngleAxis[2] };
			double q[4];
			ceres::AngleAxisToQuaternion(angle, q);


			// translation -> camera position c=-R'*t
			double t[3] = { cam->Translation[0], cam->Translation[1], cam->Translation[2] };
			double c[3];
			angle[0] = -angle[0];
			angle[1] = -angle[1];
			angle[2] = -angle[2];
			ceres::AngleAxisRotatePoint(angle, t, c);
			c[0] = -c[0];
			c[1] = -c[1];
			c[2] = -c[2];


			ofs << file << "\t"
				<< fx << " "
				<< q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " "
				<< c[0] << " " << c[1] << " " << c[2] << " "
				<< rad << " " << 0 << "\n";
		}
		ofs << endl;


		// point information
		const int nPoints = xyz.size();
		ofs << nPoints << endl;

		for (int ptID = 0; ptID < nPoints; ptID++)
		{
			double x = xyz[ptID][0];
			double y = xyz[ptID][1];
			double z = xyz[ptID][2];

			int r = 0;
			int g = 0;
			int b = 0;


			ofs << x << " " << y << " " << z << " "
				<< r << " " << g << " " << b << " ";

			for (int camID = 0; camID < nCameras; camID++)
			{
				const CameraData *cam = &camera[camID];
				vector<int>::const_iterator itr = std::find(cam->ptID.begin(), cam->ptID.end(), ptID);
				if (itr == cam->ptID.end())
					continue;

				int idx = itr - cam->ptID.begin();
				double u = cam->point2D[idx][0] - 0.5*(cam->imgWidth - 1.0);
				double v = cam->point2D[idx][1] - 0.5*(cam->imgHeight - 1.0);

				ofs << camID << " " << idx << " " << u << " " << v << " ";
			}
			ofs << endl;

		}
		ofs << "\n"
			<< "\n"
			<< "\n"
			<< "0\n"
			<< "\n"
			<< "#the last part of NVM file points to the PLY files\n"
			<< "#the first number is the number of associated PLY files\n"
			<< "#each following number gives a model - index that has PLY\n"
			<< "0"
			<< endl;

		return true;
	}
	bool saveNVM(const string path, const string outputnvmname, const vector<CameraData> &camera, const vector< vector<double> > &xyz, const vector< vector<double> > &uv)
	{
		string filename = path + "/" + outputnvmname;
		ofstream ofs(filename);
		if (ofs.fail())
		{
			cerr << "Cannot write " << filename << endl;
			return false;
		}


		// camera information
		int nCameras = camera.size();
		//int nCameras = 0;
		//for (int camID = 0; camID < camera.size(); camID++)
		//	if (camera[camID].available)
		//		nCameras++;


		ofs << "NVM_V3 \n"
			<< "\n"
			<< nCameras << endl;

		ofs << fixed << setprecision(12);

		for (int camID = 0; camID < nCameras; camID++)
		{
			const CameraData *cam = &camera[camID];
			//if (!cam->available)
			//	continue;

			string file = cam->filename;
			double fx = cam->FocalLength[0];
			double rad = -cam->Radialfirst;


			// angle axis -> quaternion
			double angle[3] = { cam->AngleAxis[0], cam->AngleAxis[1], cam->AngleAxis[2] };
			double q[4];
			ceres::AngleAxisToQuaternion(angle, q);


			// translation -> camera position c=-R'*t
			double t[3] = { cam->Translation[0], cam->Translation[1], cam->Translation[2] };
			double c[3];
			angle[0] = -angle[0];
			angle[1] = -angle[1];
			angle[2] = -angle[2];
			ceres::AngleAxisRotatePoint(angle, t, c);
			c[0] = -c[0];
			c[1] = -c[1];
			c[2] = -c[2];


			ofs << file << "\t"
				<< fx << " "
				<< q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " "
				<< c[0] << " " << c[1] << " " << c[2] << " "
				<< rad << " " << 0 << "\n";
		}
		ofs << endl;


		// point information
		const int nPoints = xyz.size();
		ofs << nPoints << endl;

		for (int ptID = 0; ptID < nPoints; ptID++)
		{
			double x = xyz[ptID][0];
			double y = xyz[ptID][1];
			double z = xyz[ptID][2];

			int r = 0;
			int g = 0;
			int b = 0;


			ofs << x << " " << y << " " << z << " "
				<< r << " " << g << " " << b << " ";

			stringstream ss;
			int count = 0;
			for (int camID = 0; camID < nCameras; camID++)
			{
				double u = uv[ptID][2 * camID];
				double v = uv[ptID][2 * camID + 1];

				if (u < 0 || v < 0)
					continue;

				const CameraData *cam = &camera[camID];

				//u -= 0.5*(cam->imgWidth  - 1.0);
				//v -= 0.5*(cam->imgHeight - 1.0);
				u -= cam->OpticalCenter[0];
				v -= cam->OpticalCenter[1];

				ss << camID << " " << ptID << " " << u << " " << v << " ";
				count++;
			}
			ofs << count << " " << ss.str() << endl;

		}
		ofs << "\n"
			<< "\n"
			<< "\n"
			<< "0\n"
			<< "\n"
			<< "#the last part of NVM file points to the PLY files\n"
			<< "#the first number is the number of associated PLY files\n"
			<< "#each following number gives a model - index that has PLY\n"
			<< "0"
			<< endl;

		return true;
	}
	void copyIntrinsic(const CameraData &src, CameraData &dst)
	{
		dst.FocalLength[0] = src.FocalLength[0];
		dst.FocalLength[1] = src.FocalLength[1];
		dst.OpticalCenter[0] = src.OpticalCenter[0];
		dst.OpticalCenter[1] = src.OpticalCenter[1];
		dst.Skew = src.Skew;

		dst.Radialfirst = src.Radialfirst;
		dst.Radialothers[0] = src.Radialothers[0];
		dst.Radialothers[1] = src.Radialothers[1];
		dst.Tangential[0] = src.Tangential[0];
		dst.Tangential[1] = src.Tangential[1];
		dst.Prism[0] = src.Prism[0];
		dst.Prism[1] = src.Prism[1];

		dst.omega = src.omega;
		dst.DistCtr[0] = src.DistCtr[0];
		dst.DistCtr[1] = src.DistCtr[1];
	}
	void copyExtrinsic(const CameraData &src, CameraData &dst)
	{
		for (int i = 0; i < 3; i++)
		{
			dst.AngleAxis[i] = src.AngleAxis[i];
			dst.Translation[i] = src.Translation[i];
		}
	}

	double MeanArray(vector<double>&data)
	{
		double mean = 0.0;
		for (int ii = 0; ii < data.size(); ii++)
			mean += data[ii];
		return mean / data.size();
	}
	double VarianceArray(vector<double>&data, double mean)
	{
		if (mean == NULL)
			mean = MeanArray(data);

		double var = 0.0;
		for (int ii = 0; ii < data.size(); ii++)
			var += pow(data[ii] - mean, 2);
		return var / (data.size() - 1);
	}
}