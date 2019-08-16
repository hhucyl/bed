#include <iostream>
#include <array>
#include <math.h>
#include <vector>
#include <list>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>

using namespace std;

double random(double ubound, double dbound, int r)//r is the intended reserved number of decimals. 
{
	if(ubound < dbound)
		cerr<<"Uper boundary must be larger than lower boundary!"<<endl;
	double ratio = pow(10.0,r);
	int u = int(ubound * ratio);
	int d = int(dbound * ratio);
	return double(d+rand()%(u-d+1))/double(ratio);
}

struct grains
{
public:
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double r = 0.0;
public://structor
	grains() = default;
	grains(double gx, double gy, double gz, double gr) :  x(gx), y(gy), z(gz), r(gr){}
public://inside function
	grains &ini(double gx, double gy, double gz, double gr)
	{
		r = gr;
		x = gx;
		y = gy;
		z = gz;
		return *this;
	}	
};

double distant(grains &g1, grains &g2)
{
	return sqrt(pow((g1.x-g2.x),2.0)+pow((g1.y-g2.y),2.0)+pow((g1.z-g2.z),2.0));
}

bool HaveContacts(grains &g, vector<grains> &bg)
{
	bool state = false;
	if (bg.size() > 0)
	{
		for(auto i = bg.begin(); i != bg.end(); ++i)
		{
			if (distant(g, *i) < (g.r+i->r))
			{	
				state = true;
				break;
			}
		}
	}
	return state;
}

void GenGrains(int total, double width, double len, double depth, double radius, vector<grains> &gs, double minh)
{
	double x = random(width-radius,radius,2);
	double y = random(len-radius,radius,2);
	double z = random(0.99*depth,0.2*depth,2);
	grains g(x,y,z,radius);
	gs.push_back(g);
	for(int i = 1; i < total; ++i)
	{	
		
		do
		{
				double x = random(width-radius,radius,2);
				double y = random(len-radius,radius,2);
				double z = random(0.99*depth,minh,2);
				g.ini(x,y,z,radius);
		}while(HaveContacts(g,gs));
		gs.push_back(g);
		cout<<i*100.0/(total*1.00)<<"%"<<endl;
	}
}

void output(FILE *os, vector<grains> &g)
{
	//fprintf(os,"x              y              r\n");
	int total = g.size();	
	fprintf(os,"%10d\n",total);
	for(auto iter : g)
		fprintf(os,"%10.5e %10.5e %10.5e %10.5e\n",iter.x,iter.y,iter.z,iter.r);
}

int main(int argc, char **argv)
{
	srand( (unsigned)time(NULL) );
	FILE *out;
	out = fopen("out3d_sphere.out","w");
	double gr =           5;
	double depth =    2000.0;
	double width =    140.0;
	double len =      50.0;
	double total =    6000.0;
	double minh = 30;
	if(argc>=2) depth = atof(argv[1]);
	vector<grains> gs;
	GenGrains(total,width,len,depth,gr,gs,minh);
	output(out,gs);
	fclose(out);
}
