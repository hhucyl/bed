#include <mechsys/flbm/Domain.h>
//#include "Readh5.h"

struct UserData
{
    double        rhomax;
    double        rhomin;
    double            dt;
    double            dx;
    Vec3_t             g;
    std::ofstream oss_ss;
	size_t Nx;
    size_t Ny;
    size_t Nz;
	char const * FileKey;      
    #ifdef USE_OCL
    cl::Buffer        bBCRho;
    cl::Buffer        bBCg;
    cl::Program       UserProgram;
    #endif
};

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    
    #ifdef USE_OCL
    if (dom.IsFirstTime)
    {
        dom.IsFirstTime = false;
        double Rho[2];
        Rho[0] = dat.rhomax;
        Rho[1] = dat.rhomin;
        dat.bBCRho      = cl::Buffer(dom.CL_Context,CL_MEM_READ_WRITE,sizeof(double)*2);
        dom.CL_Queue.enqueueWriteBuffer(dat.bBCRho,CL_TRUE,0,sizeof(double)*2,Rho);
		dom.IsFirstTime = false;
        dat.bBCg        = cl::Buffer(dom.CL_Context,CL_MEM_READ_WRITE,sizeof(cl_double3));
        cl_double3      gravity[1];
        gravity[0].s[0] = dat.g[0];
        gravity[0].s[1] = dat.g[1];
        gravity[0].s[2] = dat.g[2];
        
        dom.CL_Queue.enqueueWriteBuffer(dat.bBCg,CL_TRUE,0,sizeof(cl_double3),gravity);
        
        char* pMECHSYS_ROOT;
        pMECHSYS_ROOT = getenv ("MECHSYS_ROOT");
        if (pMECHSYS_ROOT==NULL) pMECHSYS_ROOT = getenv ("HOME");

        String pCL;
        pCL.Printf("%s/mechsys/lib/flbm/lbm.cl",pMECHSYS_ROOT);

        std::ifstream infile(pCL.CStr(),std::ifstream::in);
        std::string main_kernel_code((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
        
        std::string BC_kernel_code =
	    " void kernel ApplyGravity (global const double3 * g, global double3 * BForce, global double * Rho, global const struct lbm_aux * lbmaux) \n"
            " { \n"
                " size_t ic  = get_global_id(0); \n"
                " BForce[ic] = Rho[ic]*g[0]; \n"
            " } \n"
            /*" void kernel Left_BC (global double * RBC, global const bool * IsSolid, global double * F, global double3 * Vel, global double * Rho, global const struct lbm_aux * lbmaux) \n"
            " { \n"
                " size_t ic  = get_global_id(0); \n"
                " size_t ib  = ic*lbmaux[0].Nx; \n"
                " if (!IsSolid[ib]) \n"
                " { \n"
                    //" printf(\" %f \\n \",RBC[0]); \n"
                    " size_t iv  = ib*lbmaux[0].Nneigh; \n"
                    " F[iv+1] = 1.0/3.0*(-2*F[iv+0]-4*F[iv+10]-4*F[iv+12]-4*F[iv+14]-F[iv+2]-2*F[iv+3]-2*F[iv+4]-2*F[iv+5]-2*F[iv+6]-4*F[iv+8]+2*RBC[0]); \n"
                    " F[iv+7] = 1.0/24.0*(-2*F[iv+0]-4*F[iv+10]-4*F[iv+12]-4*F[iv+14]-4*F[iv+2] +F[iv+3]-5*F[iv+4]  +F[iv+5]-5*F[iv+6]+20*F[iv+8]+2*RBC[0]); \n"
                    " F[iv+9] = 1.0/24.0*(-2*F[iv+0]+20*F[iv+10]-4*F[iv+12]-4*F[iv+14]-4*F[iv+2]+F[iv+3]-5*F[iv+4]-5*F[iv+5]+F[iv+6]-4*F[iv+8]+2*RBC[0]); \n"
                    " F[iv+11]= 1.0/24.0*(-2*F[iv+0]-4*F[iv+10]+20*F[iv+12]-4*F[iv+14]-4*F[iv+2]-5*F[iv+3]+F[iv+4]  +F[iv+5]-5*F[iv+6]-4*F[iv+8]+2*RBC[0]); \n"
                    " F[iv+13]= 1.0/24.0*(-2*F[iv+0]-4*F[iv+10]-4 *F[iv+12]+20*F[iv+14]-4*F[iv+2]-5*F[iv+3]+  F[iv+4]-5*F[iv+5]+F[iv+6]-4*F[iv+8]+2*RBC[0]); \n"
                    " Rho   [ib] = 0.0; \n"
                    " Vel   [ib] = (double3)(0.0,0.0,0.0); \n"
                    " for(size_t k=0;k<lbmaux[0].Nneigh;k++) \n"
                    " { \n"
                        " Rho[ib] += F[iv + k]; \n"
                        " Vel[ib] += F[iv + k]*lbmaux[0].C[k]; \n"
                    " } \n"
                    " Vel[ib] *= lbmaux[0].Cs/Rho[ib]; \n"
                " } \n"
            " } \n"
            
            " void kernel Right_BC (global double * RBC, global const bool * IsSolid, global double * F, global double3 * Vel, global double * Rho, global const struct lbm_aux * lbmaux) \n"
            " { \n"
                " size_t ic  = get_global_id(0); \n"
                " size_t ib  = ic*lbmaux[0].Nx + lbmaux[0].Nx-1; \n"
                " if (!IsSolid[ib]) \n"
                " { \n"
                    //" printf(\" %f \\n \",RBC[1]); \n"
                    " size_t iv  = ib*lbmaux[0].Nneigh; \n"
                    " F[iv+2] = 1/3.0* (-2*F[iv+0]-F[iv+1]-2*(2*F[iv+11]+2*F[iv+13]+F[iv+3]+F[iv+4]+F[iv+5]+F[iv+6]+2*F[iv+7]+2*F[iv+9]-RBC[1])); \n"
                    " F[iv+8] = 1/24.0*(-2*F[iv+0] - 4*F[iv+1] - 4*F[iv+11] - 4*F[iv+13] - 5*F[iv+3] + F[iv+4] - 5*F[iv+5] + F[iv+6] +20*F[iv+7] - 4*F[iv+9] + 2*RBC[1]); \n"
                    " F[iv+10]= 1/24.0*(-2*F[iv+0] - 4*F[iv+1] - 4*F[iv+11] - 4*F[iv+13] - 5*F[iv+3] + F[iv+4] + F[iv+5] - 5*F[iv+6] - 4*F[iv+7] + 20*F[iv+9] + 2*RBC[1]) ; \n"
                    " F[iv+12]= 1/24.0*(-2*F[iv+0] - 4*F[iv+1] + 20*F[iv+11] - 4*F[iv+13] + F[iv+3] - 5*F[iv+4] - 5*F[iv+5] + F[iv+6] -  4*F[iv+7] - 4*F[iv+9] + 2*RBC[1]); \n"
                    " F[iv+14]= 1/24.0*(-2*F[iv+0] - 4*F[iv+1] - 4*F[iv+11] + 20*F[iv+13] + F[iv+3] - 5*F[iv+4] + F[iv+5] - 5*F[iv+6] -  4*F[iv+7] - 4*F[iv+9] + 2*RBC[1]); \n"
                    " Rho   [ib] = 0.0; \n"
                    " Vel   [ib] = (double3)(0.0,0.0,0.0); \n"
                    " for(size_t k=0;k<lbmaux[0].Nneigh;k++) \n"
                    " { \n"
                        " Rho[ib] += F[iv + k]; \n"
                        " Vel[ib] += F[iv + k]*lbmaux[0].C[k]; \n"
                    " } \n"
                    " Vel[ib] *= lbmaux[0].Cs/Rho[ib]; \n"
                " } \n"
            " } \n"*/
	" void kernel Top_BC (global double * RBC, global const bool * IsSolid, global double * F, global double3 * Vel, global double * Rho, global const struct lbm_aux * lbmaux) \n"
        " { \n"
	        " size_t ic  = get_global_id(0); \n"
                " size_t ib  = ic+(lbmaux[0].Nz-1)*lbmaux[0].Nx*lbmaux[0].Ny; \n"
		///" size_t ibl  = ic*lbmaux[0].Nx; \n"
		///" size_t ibr  = ic*lbmaux[0].Nx + lbmaux[0].Nx-1; \n"
		//"size_t ibx = ib%lbmaux[0].Nx;\n"
    		//"size_t iby = (ib/lbmaux[0].Nx)%lbmaux[0].Ny;\n"
    		//"size_t ibz = ib/(lbmaux[0].Nx*lbmaux[0].Ny);\n"

    		//"for (size_t k=1;k<lbmaux[0].Nneigh;k++)\n"
    		//"{\n"
        	//	"size_t inx = (size_t)((long)ibx + (long)lbmaux[0].C[k].x + (long)lbmaux[0].Nx)%lbmaux[0].Nx;\n"
        	//	"size_t iny = (size_t)((long)iby + (long)lbmaux[0].C[k].y + (long)lbmaux[0].Ny)%lbmaux[0].Ny;\n"
        	//	"size_t inz = (size_t)((long)ibz + (long)lbmaux[0].C[k].z + (long)lbmaux[0].Nz)%lbmaux[0].Nz;\n"
        	//	"size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;\n"
        	//	"F[ib*lbmaux[0].Nneigh + k] = F[in*lbmaux[0].Nneigh + k] ;\n"
    		//"}\n"	
				///" size_t ib1 = ic+(lbmaux[0].Nz-2)*lbmaux[0].Nx*lbmaux[0].Ny; \n"
				///" size_t ib2 = ic+(lbmaux[0].Nz-3)*lbmaux[0].Nx*lbmaux[0].Ny; \n"
                " if (!IsSolid[ib]) \n"
                " { \n"
                    " size_t iv = ib*lbmaux[0].Nneigh; \n"
		    		///" size_t iv1 = ib1*lbmaux[0].Nneigh; \n"
		    		///" Rho   [ib] = (4.0*Rho[ib1] - Rho[ib2])/3.0; \n"
                    ///" Vel   [ib] = (4.0*Vel[ib1] - Vel[ib2])/3.0; \n"
/*
					" bool valid = true; \n"
        			" double alpha = 1.0; \n"
        			" while (valid) \n"
        			" { \n"
            			" valid = false; \n"
		    			" for(size_t k=0; k<lbmaux[0].Nneigh; k++) \n"
		    			" { \n"
		    				" F[iv+k] = F[iv1+k] - alpha*(FeqFluid(k,Rho[ib1],Vel[ib1],lbmaux) - FeqFluid(k,Rho[ib],Vel[ib],lbmaux));\n"
							" if(F[iv+k] < -1.0e-12)"
							" { \n"
							" 	double temp = F[iv1 + k]/(FeqFluid(k,Rho[ib1],Vel[ib1],lbmaux) - FeqFluid(k,Rho[ib],Vel[ib],lbmaux));\n"
							" 	if (temp<alpha) alpha = temp;\n"
							" valid = true;"
							" } \n"
			//"F[iv+k] =  F[iv1+ lbmaux[0].Op[k]];\n"
                    	" } \n"
					 " }\n"
*/
					" F[iv+6] = F[iv+5]; \n"
					" F[iv+9] = F[iv+7] ; \n"
					" F[iv+8] = F[iv+10]; \n"
					" F[iv+13] = F[iv+11]; \n"
					" F[iv+12] = F[iv+14]; \n"
					
		    		" Rho   [ib] = 0.0; \n"
                    " Vel   [ib] = (double3)(0.0,0.0,0.0); \n"
                    " for(size_t k=0;k<lbmaux[0].Nneigh;k++) \n"
                    " { \n"
                        " Rho[ib] += F[iv + k]; \n"
                        " Vel[ib] += F[iv + k]*lbmaux[0].C[k]; \n"
                    " } \n"
                    " Vel[ib] *= lbmaux[0].Cs/Rho[ib]; \n"
				" } \n"
	" } \n"
        ;

        BC_kernel_code = main_kernel_code + BC_kernel_code;

        cl::Program::Sources sources;
        sources.push_back({BC_kernel_code.c_str(),BC_kernel_code.length()});

        dat.UserProgram = cl::Program(dom.CL_Context,sources);
        if(dat.UserProgram.build({dom.CL_Device})!=CL_SUCCESS){
            std::cout<<" Error building: "<<dat.UserProgram.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dom.CL_Device)<<"\n";
            exit(1);
        }

    }

    cl::Kernel kernel(dat.UserProgram,"ApplyGravity");
    kernel.setArg(0,dat.bBCg      );
    kernel.setArg(1,dom.bBForce[0]);
    kernel.setArg(2,dom.bRho   [0]);
    kernel.setArg(3,dom.blbmaux   );
    dom.CL_Queue.enqueueNDRangeKernel(kernel,cl::NullRange,cl::NDRange(dom.Ncells),cl::NullRange);
    dom.CL_Queue.finish();
/*
    cl::Kernel kernel(dat.UserProgram,"Left_BC");
    kernel.setArg(0,dat.bBCRho     );
    kernel.setArg(1,dom.bIsSolid[0]);
    kernel.setArg(2,dom.bF      [0]);
    kernel.setArg(3,dom.bVel    [0]);
    kernel.setArg(4,dom.bRho    [0]);
    kernel.setArg(5,dom.blbmaux    );
    dom.CL_Queue.enqueueNDRangeKernel(kernel,cl::NullRange,cl::NDRange(dom.Ndim[1]*dom.Ndim[2]),cl::NullRange);
    dom.CL_Queue.finish();

    kernel = cl::Kernel(dat.UserProgram,"Right_BC");
    kernel.setArg(0,dat.bBCRho     );
    kernel.setArg(1,dom.bIsSolid[0]);
    kernel.setArg(2,dom.bF      [0]);
    kernel.setArg(3,dom.bVel    [0]);
    kernel.setArg(4,dom.bRho    [0]);
    kernel.setArg(5,dom.blbmaux    );
    dom.CL_Queue.enqueueNDRangeKernel(kernel,cl::NullRange,cl::NDRange(dom.Ndim[1]*dom.Ndim[2]),cl::NullRange);
    dom.CL_Queue.finish();*/

    kernel = cl::Kernel(dat.UserProgram,"Top_BC");
    kernel.setArg(0,dat.bBCRho     );
    kernel.setArg(1,dom.bIsSolid[0]);
    kernel.setArg(2,dom.bF      [0]);
    kernel.setArg(3,dom.bVel    [0]);
    kernel.setArg(4,dom.bRho    [0]);
    kernel.setArg(5,dom.blbmaux    );
    dom.CL_Queue.enqueueNDRangeKernel(kernel,cl::NullRange,cl::NDRange(dom.Ndim[0]*dom.Ndim[1]),cl::NullRange);
    dom.CL_Queue.finish();


    #else // USE_OCL
    // Cells with prescribed velocity
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Ndim(1); ++i)
	for (size_t j=0; j<dom.Ndim(2); ++j)
	{
        if (dom.IsSolid[0][0][i][j]) continue;
        double * f = dom.F[0][0][i][j];
        
        f[1] = 1.0/3.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-f[2]-2*f[3]-2*f[4]-2*f[5]-2*f[6]-4*f[8]+2*dat.rhomax);
        f[7] = 1.0/24.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-4*f[2] +f[3]-5*f[4]  +f[5]-5*f[6]+20*f[8]+2*dat.rhomax);
        f[9] = 1.0/24.0*(-2*f[0]+20*f[10]-4*f[12]-4*f[14]-4*f[2]+f[3]-5*f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);
        f[11]= 1.0/24.0*(-2*f[0]-4*f[10]+20*f[12]-4*f[14]-4*f[2]-5*f[3]+f[4]  +f[5]-5*f[6]-4*f[8]+2*dat.rhomax);
        f[13]= 1.0/24.0*(-2*f[0]-4*f[10]-4 *f[12]+20*f[14]-4*f[2]-5*f[3]+  f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);

        dom.Vel[0][0][i][j] = OrthoSys::O;
        dom.Rho[0][0][i][j] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][0][i][j] +=  dom.F[0][0][i][j][k];
            dom.Vel[0][0][i][j] +=  dom.F[0][0][i][j][k]*dom.C[k];
        }
        dom.Vel[0][0][i][j] *= dom.Cs/dom.Rho[0][0][i][j];
	}

	// Cells with prescribed density
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Ndim(1); ++i)
	for (size_t j=0; j<dom.Ndim(2); ++j)
	{
        if (dom.IsSolid[0][dom.Ndim(0)-1][i][j]) continue;
        double * f = dom.F[0][dom.Ndim(0)-1][i][j];

        f[2] = 1/3.0* (-2*f[0]-f[1]-2*(2*f[11]+2*f[13]+f[3]+f[4]+f[5]+f[6]+2*f[7]+2*f[9]-dat.rhomin));
        f[8] = 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] - 5*f[5] + f[6] +20*f[7] - 4*f[9] + 2*dat.rhomin);
        f[10]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] + f[5] - 5*f[6] - 4*f[7] + 20*f[9] + 2*dat.rhomin) ;
        f[12]= 1/24.0*(-2*f[0] - 4*f[1] + 20*f[11] - 4*f[13] + f[3] - 5*f[4] - 5*f[5] + f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        f[14]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] + 20*f[13] + f[3] - 5*f[4] + f[5] - 5*f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        
        dom.Vel[0][dom.Ndim(0)-1][i][j] = OrthoSys::O;
        dom.Rho[0][dom.Ndim(0)-1][i][j] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][dom.Ndim(0)-1][i][j] +=  dom.F[0][dom.Ndim(0)-1][i][j][k];
            dom.Vel[0][dom.Ndim(0)-1][i][j] +=  dom.F[0][dom.Ndim(0)-1][i][j][k]*dom.C[k];
        }
        dom.Vel[0][dom.Ndim(0)-1][i][j] *= dom.Cs/dom.Rho[0][dom.Ndim(0)-1][i][j];
	}
    #endif // USE_OCL
}

void Report (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
	String fn;
	fn.Printf    ("%s_%04d", dat.FileKey, dom.idx_out);
    fn.append(".h5");

    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	size_t  Nx = dom.Ndim(0);
    size_t  Ny = dom.Ndim(1);
    size_t  Nz = dom.Ndim(2);	
	size_t  Nneigh = dom.Nneigh;
    double *Ff   = new double[Nneigh*Nx*Ny*Nz];
	size_t i = 0;
	for (size_t m=0;m<Nz;m++)
    for (size_t l=0;l<Ny;l++)
    for (size_t n=0;n<Nx;n++)
    {
        for (size_t k=0; k<Nneigh; k++)
        {
			Ff[Nneigh*i + k] = (double) dom.F[0][n][l][m][k];
        }
        i++;
    }
        
    //Writing data to h5 file
    hsize_t dims[1];
    dims[0] = Nneigh*Nx*Ny*Nz;
    String dsname;
    dsname.Printf("F_%d",0);
    H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ff);

    delete [] Ff;
	H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

}




void InitFromH5(FLBM::Domain &dom, char const * FileKey1, char const * FileKey2)
{
	hid_t file_id1;
	hid_t file_id2;
	file_id1 = H5Fopen(FileKey1,H5F_ACC_RDONLY,H5P_DEFAULT);
	file_id2 = H5Fopen(FileKey2,H5F_ACC_RDONLY,H5P_DEFAULT);
	double *Ff = new double[dom.Nneigh*dom.Ndim(0)*dom.Ndim(1)*dom.Ndim(2)];
	double *vel = new double[3*dom.Ndim(0)*dom.Ndim(1)*dom.Ndim(2)];
	double *rho = new double[dom.Ndim(0)*dom.Ndim(1)*dom.Ndim(2)];
	H5LTread_dataset_double(file_id1,"/F_0",Ff);
	H5LTread_dataset_double(file_id2,"/Velocity_0",vel);
	H5LTread_dataset_double(file_id2,"/Density_0",rho);
	size_t nn=0;
	for (size_t iz=0; iz<dom.Ndim(2); ++iz)
    for (size_t iy=0; iy<dom.Ndim(1); ++iy)
    for (size_t ix=0; ix<dom.Ndim(0); ++ix)
    {
		double * f  = dom.F[0][ix][iy][iz];
		for(size_t k=0; k<dom.Nneigh; k++)
		{
			
			f[k] = Ff[dom.Nneigh*nn + k];
		}
		dom.Vel[0][ix][iy][iz] = vel[3*nn],vel[3*nn+1],vel[3*nn+2];
        dom.Rho[0][ix][iy][iz] = rho[nn];
		nn++;
	}

	delete[] Ff;
	H5Fflush(file_id1,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id1);

	H5Fflush(file_id2,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id2);
	
}

void InitFromMsys(FLBM::Domain &dom, double rho, Vec3_t &vv)
{
	size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
	// Data data;
	// const char *fn = "/home/user/hyporheic_stimulation/t3d_top_experiment/feh/Feh_300_0097.h5";
	// data.Readh5(fn);	
    for (size_t ix=0; ix<nx; ++ix)
    for (size_t iy=0; iy<ny; ++iy)
    for (size_t iz=0; iz<nz; ++iz)
    {
		// v0 = data.Velocity[3*(i-1)], data.Velocity[3*(i-1)+1], data.Velocity[3*(i-1)+2];
		// rho0 = data.Density[i];
		// if(Dom.Lat[0].Cells[i]->IsSolid)
		// {
		// 	Dom.Lat[0].Cells[i]->Initialize(1.0,vv);
		// }else{
		// 	if(rho0 < 0.7)
		// 	{
		// 		Dom.Lat[0].Cells[i]->Initialize(1.0 , v0);
		// 	}else{
		// 		Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
		// 	}
		// }	
        Vec3_t idx(ix,iy,iz);
        dom.Initialize(0,idx,rho,vv);
    }
}
void AddSphereQ(FLBM::Domain &dom, Vec3_t &pos, double R)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    double px = pos(0);
    double py = pos(1);
    double pz = pos(2);
    int xstart = std::max(px-R-2,0.0);
    int ystart = std::max(py-R-2,0.0);
    int zstart = std::max(pz-R-2,0.0);
    int xend = std::min(px+R+2,(double) nx);
    int yend = std::min(py+R+2,(double) ny);
    int zend = std::min(pz+R+2,(double) nz);

    //#ifdef USE_OMP
    //#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    //#endif
    for(int ix= xstart; ix<xend; ++ix)
    for(int iy= ystart; iy<yend; ++iy)
    for(int iz= zstart; iz<zend; ++iz)  
    {
        double x = (double) ix;
        double y = (double) iy;
        double z = (double) iz;
        Vec3_t Cell(x,y,z);
        double L = dot(Cell-pos,Cell-pos);
        if(L<R*R) 
        {
            dom.IsSolid[0][ix][iy][iz] = true;
            continue;
        }
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            size_t nix = (size_t)((int)ix + (int)dom.C[k](0) + (int)nx)%nx;
            size_t niy = (size_t)((int)iy + (int)dom.C[k](1) + (int)ny)%ny;
            size_t niz = (size_t)((int)iz + (int)dom.C[k](2) + (int)nz)%nz;
            x = (double) nix;
            y = (double) niy;
            z = (double) niz;
            Vec3_t Cell1(x,y,z);
            double L1 = dot(Cell1-pos,Cell1-pos);
            if(L1>R*R) continue;
            // Gamma[ix][iy][iz] = 1.0;                
            
            double K = dot(Cell-pos,Cell1-pos);
            double A = (L + L1- 2.0*K);
            double B = (2.0*K - 2.0*L);
            double C = L - R*R;
            double delta = B*B - 4.0*A*C;
            if(delta<0)
            {
                continue;
            }
            double q1 = (-B + std::sqrt(delta))/(2.0*A);
            double q2 = (-B - std::sqrt(delta))/(2.0*A);
            bool flag1 = q1>=0 && q1-1 <1e-6;
            bool flag2 = q2>=0 && q2-1 <1e-6;
            if(flag1)
            {
                dom.IsSolid[0][nix][niy][niz] = true;
            }else{
                if(flag2)
                {
                    dom.IsSolid[0][nix][niy][niz] = true;
                }

            }       
        }
        
    }
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1;
    size_t h = 5;
	size_t N = 5;
    double Re = 39000;	
    double ua = 0.24;
    double H = 0.160;
    double ttt = 1e-4;//t/tlbm
    double Scc = 0.17;
    
    
    if (argc>=2) h = atoi(argv[1]);
    if (argc>=3) Scc   =atof(argv[2]);
    if (argc>=4) ttt   =atof(argv[3]);
    if (argc>=5) Re   =atof(argv[4]);
    double R = (double) h;
    size_t nx = 2*N*h;
    size_t ny = 10*h+2;
    size_t nz = 12*h+2+std::ceil(H/0.038*2*h);//0.12/0.038 = 3.1579
    double lll = 0.038/((double) 2*h);//L/Llbm
    double dx = 1.0;
    double dt = 1.0;
    double nu = H*ua/Re*ttt/(lll*lll);

    Vec3_t pos(0,0,0);

    FLBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    
    UserData dat;
    Dom.UserData = &dat;
    Dom.Sc = Scc;

    dat.rhomin  = 1.0;
    dat.rhomax  = 1.000001;
    dat.dx      = dx;
    dat.dt      = dt;
    dat.g           = 1.5e-4,0.0,0.0;
	dat.Nx = nx;
    dat.Ny = ny;
    dat.Nz = nz;
	dat.FileKey = "perbed_F"; 
    //Assigning solid boundaries at bottom
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        Dom.IsSolid[0][i][j][0]    = true;
		Dom.IsSolid[0][i][j][1]    = true;
        //Dom.IsSolid[0][i][j][nz-1] = true;
    }
	for(size_t ix=0; ix<nx; ++ix)
	for(size_t iz=0; iz<nz; ++iz)
	{
		Dom.IsSolid[0][ix][0][iz] = true;
		Dom.IsSolid[0][ix][ny-1][iz] = true;
	}

    for(size_t ipy=0; ipy<5; ipy++)
    for(size_t ipz=0; ipz<6; ipz++)
    for(size_t ipx=0; ipx<N; ipx++)
    {
        pos = (2*ipx+1)*R,(2*ipy+1)*R + 1,(2*ipz+1)*R + 2;
        AddSphereQ(Dom,pos,R);
    }
	Vec3_t vv(0.0,0.0,0.0);
	InitFromMsys(Dom,1.0,vv);
	//InitFromH5(Dom,"perbed_F_0001.h5","perbed_0001.h5");
    //Dom.G [0]    = -200.0;
    //Dom.Gs[0]    = -400.0;
    //Dom.Gmix     =  0.001;  

    Dom.Solve(500*dt,1.0e2*dt,Setup,Report,"perbed",true,Nproc);
    //Dom.Solve(1.0e2,1,Setup,NULL,"single4",true,Nproc);
    //Dom.Solve(1.0,80.0,NULL,NULL,"single",true,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH
