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
	double Tlen = 50.0;
	double dlen = 100.0;
	double depth = 2100.0;
	double width = 140.0;

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
	Vec3_t v(random(0.1,-0.1,6),random(0.1,-0.1,6),-0);
	Vec3_t w(0.1,0.1,0.1);
	std::fstream ifile("out3d_sphere.out",std::ios::in);

	if(!ifile.fail())
	{
		//while(!ifile.eof())
		int total;
		ifile>>total;
		for(int i =0; i< total; ++i)
		{
		double x;
    		double y;
    		double z;
    		double radius;
			//ifstream ParData ("particles45.txt");
    		ifile>>x>>y>>z>>radius;
		Vec3_t xx(x,y,z);
		d.AddSphere(-6, xx, radius, 0.1);
		d.Particles.Last()->v = v;
		d.Particles.Last()->w = w;
		d.Particles.Last()->Ff = (d.Particles.Last()->Props.m)*Vec3_t(0.0,0.0,-0.5);
		d.Particles.Last()->Props.Gn = -0.2;
		}
	}		
	

	//set properties    
	Dict B;
    	for(int i=1; i < 6; ++i)
	{
		B.Set(-i,"Gt mu Kn Kt",0.0,0.4,1.0e8,5.0e7);
		
	}
	d.SetProps(B);
	// solve
        d.CamPos = 0.0,30.0,0.0;
        d.Solve(/*tf*/200.0, 1.0e-4, /*dtOut*/4, NULL, NULL, "td_more", 2, 12);
v=0.0,0.0,-10.0;
	w=0.0,0.0,0.0;
	//pressure plane
	xx=0.5*width,0.5*Tlen,1500;
	d.AddPlane(-7,xx,10.,width,Tlen,10.,0,&OrthoSys::e0);
	d.GetParticle(-7)->FixVeloc();	
	d.GetParticle(-7)->v = Vec3_t(0.0,0.0,-10.0);
	d.GetParticle(-7)->w = Vec3_t(0.0,0.0, 0.0);
	d.GetParticle(-7)->vzf = true;
	B.Set(-7,"Gt mu Kn Kt",0.0,0.4,1.0e8,5.0e7);
	d.SetProps(B);
	d.CamPos = 0.0,30.0,0.0;
        d.Solve(/*tf*/210.0, 1.0e-4, /*dtOut*/1, NULL, NULL, "td_more2", 2, 12);
}MECHSYS_CATCH


