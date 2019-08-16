// MechSys
#include <mechsys/dem/domain.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>

double random(double ubound, double dbound, int r)//r is the intended reserved number of decimals. 
{
	if(ubound < dbound)
		std::cerr<<"Uper boundary must be larger than lower boundary!"<<std::endl;
	double ratio = pow(10.0,r);
	int u = int(ubound * ratio);
	int d = int(dbound * ratio);
	return double(d+rand()%(u-d+1))/double(ratio);
}
int main(int argc, char **argv) try
{
	//random number
	srand( (unsigned)time(NULL) );
	double Tlen = 150.0;
	double dlen = 100.0;
	double depth = 2100.0;
	double width = 300.0;

	DEM::Domain d;
	//xz1
	Vec3_t xx(0.5*width,0.0,depth);
	d.AddPlane(-1,xx,1,width,2*depth,1,-M_PI/2,&OrthoSys::e0);
	//yz1	
	xx=0,0.5*Tlen,depth;
	d.AddPlane(-2,xx,1,2.01*depth,Tlen,1,-M_PI/2,&OrthoSys::e1);
	//xz2	
	xx=0.5*width,Tlen,depth;
	d.AddPlane(-3,xx,1,width,2*depth,1,-M_PI/2,&OrthoSys::e0);
	//yz2	
	xx=width,0.5*Tlen,depth;	
	d.AddPlane(-4,xx,1,2.01*depth,Tlen,1,-M_PI/2,&OrthoSys::e1);
	//xy
	xx=0.5*width,0.5*Tlen,0.0;
	d.AddPlane(-5,xx,1,width,Tlen,1,0,&OrthoSys::e0);
	
	for(int i=1; i<6; i++)
	{	
		d.GetParticle(-i)->FixVeloc();
	}
	//generate the grains
	size_t box_num = 5;
	std::fstream ifile("out3d_sphere.out",std::ios::in);
	int total;
	ifile>>total;
	char const *FileKey;
	FileKey = "td_more_0049.h5";
	hid_t file_id;
	file_id = H5Fopen(FileKey,H5F_ACC_RDONLY,H5P_DEFAULT);
	double *Pos;
	Pos = new double[3*(box_num+total)];
	H5LTread_dataset_double(file_id,"/Position",Pos);
	double *R;
	R = new double[box_num+total];
	H5LTread_dataset_double(file_id,"/Radius",R);
	
	
	for(int i =0; i< total; ++i)
	{
		size_t ix = 3*box_num+3*i;
		size_t iy = 3*box_num+3*i+1;
		size_t iz = 3*box_num+3*i+2;
		Vec3_t xx(Pos[ix],Pos[iy],Pos[iz]);
		std::cout<<Pos[ix]<<" "<<Pos[iy]<<" "<<Pos[iz]<<" "<<R[box_num+i]<<std::endl;
		if(Pos[iz]>300.0) continue;
		d.AddSphere(-6, xx, R[box_num+i], 0.1);
		d.Particles.Last()->v = Vec3_t(0.0,0.0,0.0);
		d.Particles.Last()->w = Vec3_t(0.0,0.0,0.0);
		d.Particles.Last()->Ff = (d.Particles.Last()->Props.m)*Vec3_t(0.0,0.0,-0.5);
		d.Particles.Last()->Props.Gn = -0.2;
		}
		
	

	//set properties    
	Dict B;
    	for(int i=1; i < 6; ++i)
	{
		B.Set(-i,"Gt mu Kn Kt",0.0,0.4,1.0e10,5.0e7);
		
	}
	d.SetProps(B);
	
	//pressure plane
	xx=0.5*width,0.5*Tlen,312;
	d.AddPlane(-7,xx,10.,width,Tlen,10.,0,&OrthoSys::e0);
	d.GetParticle(-7)->FixVeloc();	
	d.GetParticle(-7)->v = Vec3_t(0.0,0.0,-1.0);
	d.GetParticle(-7)->w = Vec3_t(0.0,0.0, 0.0);
	d.GetParticle(-7)->vzf = true;
	d.GetParticle(-5)->FixVeloc();	
	d.GetParticle(-5)->v = Vec3_t(0.0,0.0,1.0);
	d.GetParticle(-5)->w = Vec3_t(0.0,0.0, 0.0);
	d.GetParticle(-5)->vzf = true;
	B.Set(-7,"Gt mu Kn Kt",0.0,0.4,1.0e10,5.0e7);
	d.SetProps(B);
	d.CamPos = 0.0,30.0,0.0;
        d.Solve(/*tf*/500.0, 1.0e-4, /*dtOut*/1, NULL, NULL, "td_more2", 2, 1);
}MECHSYS_CATCH


