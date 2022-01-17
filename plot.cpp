#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "gnuplot-iostream.h"
using namespace std;

int main()
{
	class field_data
	{
	public:
		double x;
		double y;
		double z;
		double u;
		double v;
		double w;
		double p;
	};
	
	ifstream file("/home/john/results_fyp/pipeflow_same_with_relaxation/no_head.dat");
	int plot = 0;
	
	
	if(file.fail())
	{
		cout<<"Error opening file"<<endl;
		exit(1);
	}
	
	
	vector<field_data> f;

	 int count = 0;
	 
	while(!file.eof())
	{
		field_data temp;
		file>>temp.x>>temp.y>>temp.z>>temp.u>>temp.v>>temp.w>>temp.p;
		f.push_back(temp);
	}
	double sum = 0.0;
	for (unsigned int i = 0; i<f.size(); i++)
	{
		if(f[i].x== 0.00206667 and f[i].y == 0.00206667 and f[i].z> 0.44)
		{
			sum += (f[i].p - f[i-1].p)/(f[i-1].z - f[i].z);
			count++;
		}
	}
	cout<<" p-grad avg"<<sum/count<<endl;
//	cout<<sum<<endl;
/*
	double max_x = 0.0;
	double max_y = 0.0;
	double max_w = 0.0;
	
	for (unsigned int i = 0; i<f.size(); i++)
	{
		if(f[i].z== 0.03)
		{
			if(f[i].w >= max_w)
			{
				max_w = f[i].w;
				max_x  = f[i].x;
				max_y  = f[i].y;
			}
		}
	}	

	for (unsigned int i = 0; i<f.size(); i++)
	{
		if(f[i].x == max_x and f[i].y == )
		{
			cout<<f[i].w<<endl;
		}
	}	
	
	cout<<count<<endl;
	*/
	
	vector<std::pair<double,double>> data0;
	vector<std::pair<double,double>> data1;
	vector<std::pair<double,double>> data2;
	vector<std::pair<double,double>> data3;
	vector<std::pair<double, double>> data4;
	vector<std::pair<double,double>> data5;
	vector<std::pair<double, double>> data6;
	vector<std::pair<double, double>> data7;
	vector<std::pair<double, double>> data8;
	vector<std::pair<double, double>> data9;
	vector<std::pair<double, double>> data10;
	vector<std::pair<double, double>> datazu;
	
	for (unsigned int i = 0; i<f.size(); i++)
	{
		

		if(f[i].z== 0.02298013 and f[i].y == 0.0001)
		{
					data0.push_back(make_pair(f[i].x, f[i].w));
		}	

		if(f[i].z== 0.02456954 and f[i].y == 0.0001)
		{
					data1.push_back(make_pair(f[i].x, f[i].w));
		}	
		
		if(f[i].z== 0.02562914 and f[i].y == 0.0001)
		{
					data2.push_back(make_pair(f[i].x, f[i].w));
		}	
		
		if(f[i].z== 0.02655629 and f[i].y == 0.0001)
		{
					data3.push_back(make_pair(f[i].x, f[i].w));
		}	
		
		if(f[i].z== 0.02761589 and f[i].y == 0.0001)
		{
					data4.push_back(make_pair(f[i].x, f[i].w));
		}	
		
		if(f[i].z== 0.03000000 and f[i].y == 0.0001)
		{
					data5.push_back(make_pair(f[i].x, f[i].w));
		}
	
		if(f[i].z== 0.03701987 and f[i].y == 0.0001)
		{
					data6.push_back(make_pair(f[i].x, f[i].w));
		}	
		
		if(f[i].z == 0.04 and f[i].y == 0.0001)
		{
					data7.push_back(make_pair(f[i].x, f[i].w));
		}
		
		if(f[i].z == 0.04504983 and f[i].y == 0.0001)
		{
					data8.push_back(make_pair(f[i].x, f[i].w));
		}
		
		if(f[i].z == 0.05009967 and f[i].y == 0.0001)
		{
					data9.push_back(make_pair(f[i].x, f[i].w));
		}
		
		if(f[i].z == 0.05807309 and f[i].y == 0.0001)
		{
					data10.push_back(make_pair(f[i].x, f[i].w));
		}
		
		if(f[i].x == 0.0001 and f[i].y == 0.0001)
		{
					datazu.push_back(make_pair(f[i].z, f[i].u));
		}
	}	
	
if(plot==1)
{

  Gnuplot gp;
  gp << "set terminal png enhanced size 1280, 1024\n";
  gp << "set output '/home/john/results_fyp/microchannel/0.2mmx0.2mmx20mm/re=110_rew=5.02_rps=125.5/Re=110_Rew=5.02_2.png'\n";
  gp << "set xrange [0.0:0.0002]\n";
  gp << "set title 'Rew = 5.02, Re = 110, rps = 125.5, Win = 0.55m/s'\n";
  gp << "plot '-' with lines title 'Z = 0.022', '-' with points title 'Z = 0.024', '-' with points title 'Z = 0.025','-' with lines title 'Z = 0.026', '-' with points title 'Z = 0.027', '-' with points title 'Z = 0.030', '-' with points title 'Z = 0.037', '-' with points title 'Z = 0.040'\n";
  gp.send1d(data0);
  gp.send1d(data1);
  gp.send1d(data2);
  gp.send1d(data3);
  gp.send1d(data4);
  gp.send1d(data5);
  gp.send1d(data6);
  gp.send1d(data7);
  /*
  gp.send1d(data8);
  gp.send1d(data9);
  gp.send1d(data10);
  */
 }
 /*
  Gnuplot gp1;
  gp1 << "set terminal png enhanced size 1280, 1024\n";
  gp1 << "set output '/home/john/results_fyp/microchannel/0.2mmx0.2mmx40mm/re=110_rew=5.02_rps=125.5/Re=110_Rew=5.02_UvsZ.png'\n";
  gp1 << "set xrange [0.02:0.06]\n";
  gp1 << "set title 'Rew = 5.02, Re = 110, rps = 125.5, Win = 0.55m/s'\n";
  gp1 << "plot '-' with lines title 'X=0.0001  Y = 0.0001'\n";
  gp1.send1d(datazu);
 */
  return 0;
}
