/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/


//STD
#include<iostream>

// MechSys
#include <mechsys/flbm/Domain.h>
//#include <mechsys/lbm/Dompargen.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>
#include <fstream>
#include <math.h>
//my
// #include "Readh5.h"
using namespace std;
struct UserData
{
    size_t Nx;
    size_t Ny;
    size_t Nz;
	double        rhomax;
    double        rhomin;
    double            dt;
    double            dx;
    Vec3_t             g;
    std::ofstream oss_ss;
	char const * FileKey;      
   
};

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

void Setup (FLBM::Domain & dom, void * UD)
{
	UserData & dat = (*static_cast<UserData *>(UD));

    	// TOP BOUNDARY CONDITION

 	#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t ix=0; ix<dat.Nx; ++ix)
    for (size_t iy=0; iy<dat.Ny; ++iy)
	{	
		double * f  = dom.F[0][ix][iy][dat.Nz-1];
        f[6] = f[5]; 
        f[13] = f[11]; 
        f[12] = f[14]; 
        f[17] = f[15];
        f[16] = f[18];
        dom.Vel[0][ix][iy][dat.Nz-1] = 0.0, 0.0, 0.0;
        dom.Rho[0][ix][iy][dat.Nz-1] = 0.0;
        for (size_t k=0; k<dom.Nneigh; ++k)
        {
            dom.Rho[0][ix][iy][dat.Nz-1] += dom.F[0][ix][iy][dat.Nz-1][k];
            dom.Vel[0][ix][iy][dat.Nz-1] += dom.F[0][ix][iy][dat.Nz-1][k]*dom.C[k];
        }
        dom.Vel[0][ix][iy][dat.Nz-1] /= dom.Rho[0][ix][iy][dat.Nz-1];

	}

    //body force
    
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t ix=0; ix<dat.Nx; ++ix)
    for (size_t iy=0; iy<dat.Ny; ++iy)
    for (size_t iz=0; iz<dat.Nz; ++iz)
    {
        dom.BForce[0][ix][iy][iz] = dom.Rho[0][ix][iy][iz]*dat.g;

    }
	
//THE ACCERLARATION
   /*#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
   for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c   = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.g;
    }*/
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
    size_t Nproc = 12;
    size_t h = 10;
	size_t N = 12;
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

    FLBM::Domain Dom(D3Q19, nu, iVec3_t(nx,ny,nz), dx, dt);
    
    UserData dat;
    Dom.UserData = &dat;
    Dom.Sc = Scc;

    dat.rhomin  = 1.0;
    dat.rhomax  = 1.000001;
    dat.dx      = dx;
    dat.dt      = dt;
    dat.g           = 1.5e-8,0.0,0.0;
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

