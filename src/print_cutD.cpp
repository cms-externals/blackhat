#include <FeynDiagram/fd.h>
#include <fstream>
#include <valarray>
#include <algorithm>
#include "print_cutD.h"
#include <stdio.h>
#include <stdlib.h>
#include "partitions.h"
#include "cut_part.h"
#include "rec_tree.h"
#include "cut_part_normal.h"
#include "ratext/ratext_part_normal.h"
#include <stdio.h>
#include <stdlib.h>

namespace BH {



void print_box_graph(boxD* cd,const char* path);
void print_triangle_graph(triangleD* cd,const char* path);
void print_bubble_graph(bubbleD* cd,const char* path);

long cutD_external_code(const cutD& cd){
	long base=10;
	long factor=1;
	long res=0;
	for (int cor=cd.nc();cor>=1;cor--){
		for (int i=cd.c(cor).size()-1;i>=0;i--){
			res+=factor*cd.c(cor)[i].ind();
			factor*=base;
		}
	}
	return res;
}

long raw_part_external_code(const raw_part& cd){
	long base=10;
	long factor=1;
	long res=0;
	for (int cor=cd.nc();cor>=1;cor--){
		for (int i=cd.indices(cor).size()-1;i>=0;i--){
			res+=factor*cd.indices(cor)[i];
			factor*=base;
		}
	}
	return res;
}


void add_line(FeynDiagram& fd,vector<line_raw*>& lr,const particle_ID& p,const xy& start , const xy& end ){
	if (p.is_a(gluon)){
		  lr.push_back(new line_spring(fd,start,end));
	} else
	if (p.is_a(quark) || p.is_a(gluino)){
		if (((p.flavor()+100)%2) == 1 ){
			if (p.is_anti_particle()) {
				lr.push_back(new line_plain(fd,end,start));
			}
			else {
				lr.push_back(new line_plain(fd,start,end));
			}
		}
		else {
			if (p.is_anti_particle()) {
				lr.push_back(new line_doubleplain(fd,end,start));
			}
			else {
				lr.push_back(new line_doubleplain(fd,start,end));
			}
		}
		if (p.is_a(gluino)){
			lr.back()->dashon.settrue();
			lr.back()->thickness.scale(2.);
			lr.back()->dsratio.set(4);
		}
	} else
	if (p.is_a(lepton)){
		lr.push_back(new line_plain(fd,start,end)); lr.back()->dashon.settrue(); lr.back()->arrowon.setfalse();
	} else
	if (p.is_a(scalar)){
		lr.push_back(new line_doubleplain(fd,start,end));
		lr.back()->dashon.settrue();
		lr.back()->arrowon.setfalse();
		lr.back()->thickness.scale(2.);
	} else
	if (p.is_a(photon)){
		lr.push_back(new line_wiggle(fd,start,end));
	} else
	if (p.is_a(quark_massive)){
			  if (p.is_anti_particle()) {
				  lr.push_back(new line_doubleplain(fd,end,start));
			  }
			  else {
				  lr.push_back(new line_doubleplain(fd,start,end));
			  }
			  lr.back()->thickness.scale(2.);
	} else
	if (p.is_a(gluon_massive)||p.is_a(scalar_massive)){
		  lr.push_back(new line_plain(fd,start,end));
		  lr.back()->dashon.settrue(); lr.back()->arrowon.setfalse();lr.back()->thickness.scale(2.);
	} else {
		_WARNING3("Unknown type in add_line for type: ",p.type()->name(),".");
		if (((p.flavor()+100)%2) == 1 ){
			  if (p.is_anti_particle()) {
				  lr.push_back(new line_doubleplain(fd,end,start));
			  }
			  else {
				  lr.push_back(new line_doubleplain(fd,start,end));
			  }

		  }
			lr.back()->thickness.scale(3.);
	}

}


void print_box(raw_box* thebox,const char* path,int index=0){
	double factor=2;
	double f_long=2.8;
	double f_short=1.1;
	page pg;
	FeynDiagram fd(pg);

	xy c[4]={ xy(-factor,factor), xy(factor,factor), xy(factor,-factor), xy(-factor,-factor)};
vertex_dot v1(fd,c[0]), v2(fd,c[1]), v3(fd,c[2]), v4(fd,c[3]);

	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][2]={"x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x"};
	  char helicities[8][2]={"1","2","3","4","5","6","7","8"};
	  xy hel_pos[8]={ xy(-f_long,-f_short),xy(-f_long,f_short), xy(-f_short,f_long), xy(f_short,f_long), xy(f_long,f_short),xy(f_long,-f_short), xy(f_short,-f_long), xy(-f_short,-f_long)};
	  for (int i=1;i<=4;i++){
		  add_line(fd,lr,lp,c[i-1],c[(i+2)%4]);

	}
	  int nbr_ext;



	int totallabel=0;

	for (int l=1;l<=4;l++){
		nbr_ext=thebox->indices(l).size();
		vector<xy> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double angle=5*M_PI/4.-(l*M_PI/2.)+(nbr_ext-2*i+1.)/double(nbr_ext) * M_PI/4.;
			e4[i-1]=c[l-1]+xy(2*cos(angle)*factor,2*sin(angle)*factor);
			index4[i-1]=c[l-1]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);

			add_line(fd,lr,thebox->c(l)[i-1],c[l-1],e4[i-1]);

			sprintf(temps[totallabel],"%d",thebox->indices(l)[i-1]);
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}

	  if (index!=0){
		  std::ostringstream o;
		  o << index;
		  char temp[]="xxx";
			sprintf(temp,"%d",index);

		t.push_back(new text(fd,temp,xy(0,0),1.,0.5));
	  }


	std::ofstream outfile;

	outfile.open(path);
	  pg.output(outfile);
	  outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }
  }


void print_pentagon(pentagonD* thepent,const char* path,int index=0){
	double factor=2;
	double f_long=2.8;
	double f_short=1.1;
	  page pg;
	  FeynDiagram fd(pg);

	  double pos_x[5];
	  double pos_y[5];

	  for (int i =0;i<5;i++){
		  double co=cos(6.2830/4.-i*6.2830/5.);
		  double si=sin(6.2830/4.-i*6.2830/5.);
		  pos_x[i]=factor*co*2.0;
		  pos_y[i]=factor*si*2.0;
	  }

	  xy c[5]={
			  xy(pos_x[0],pos_y[0]),
			  xy(pos_x[1],pos_y[1]),
			  xy(pos_x[2],pos_y[2]),
			  xy(pos_x[3],pos_y[3]),
			  xy(pos_x[4],pos_y[4])
	  };

	  vertex_dot v1(fd,c[0]), v2(fd,c[1]), v3(fd,c[2]), v4(fd,c[3]), v5(fd,c[4]);

	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[25][3]={"xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx"};
	  char helicities[10][2]={"1","2","3","4","5","6","7","8","9","0"};

	  double hel_pos_x[10];
	  double hel_pos_y[10];

	  double cc[]={-1.,1.};
	  for (int i =0;i<5;i++){
		  for (int j=0;j<2;j++){
			  double co=cos(6.2830/4.-i*6.2830/5.-cc[j]*6.2830/20.);
			  double si=sin(6.2830/4.-i*6.2830/5.-cc[j]*6.2830/20.);
			  hel_pos_x[2*i+j]=factor*2.2*co;
			  hel_pos_y[2*i+j]=factor*2.2*si;
		  }

	  }

	  xy hel_pos[10]={
			  xy(hel_pos_x[0],hel_pos_y[0]),
			  xy(hel_pos_x[1],hel_pos_y[1]),
			  xy(hel_pos_x[2],hel_pos_y[2]),
			  xy(hel_pos_x[3],hel_pos_y[3]),
			  xy(hel_pos_x[4],hel_pos_y[4]),
			  xy(hel_pos_x[5],hel_pos_y[5]),
			  xy(hel_pos_x[6],hel_pos_y[6]),
			  xy(hel_pos_x[7],hel_pos_y[7]),
			  xy(hel_pos_x[8],hel_pos_y[8]),
			  xy(hel_pos_x[9],hel_pos_y[9])
	  };

	  for (int i=1;i<=5;i++){
		  add_line(fd,lr,thepent->l(i),c[i-1],c[(i+3)%5]);

		sprintf(helicities[2*i-2],"%s",(thepent->l(i).helicity()==-1)?"+":"-");
		sprintf(helicities[2*i-1],"%s",(thepent->l(i).helicity()==-1)?"-":"+");
		t.push_back(new text(fd,helicities[2*i-2],hel_pos[(i*2+7)%10],0.5,0.5));
		t.push_back(new text(fd,helicities[2*i-1],hel_pos[2*i-2],0.5,0.5));

	}
	  int nbr_ext;



	int totallabel=0;

	for (int l=1;l<=5;l++){
		nbr_ext=thepent->c(l).size();
		vector<xy> e5(nbr_ext);
		vector<xy> index5(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double angle=9.0*M_PI/10.-(l*2.*M_PI/5.)  -(nbr_ext-2*i+1.)/double(nbr_ext) * M_PI/5.;
			e5[i-1]=    c[l-1]+xy(2  *cos(angle)*factor,2  *sin(angle)*factor);
			index5[i-1]=c[l-1]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);

			add_line(fd,lr,thepent->c(l)[i-1],c[l-1],e5[i-1]);

			sprintf(temps[totallabel],"%d%s",thepent->c(l)[i-1].ind(),(thepent->c(l)[i-1].helicity()==-1)?"-":"+");
			t.push_back(new text(fd,temps[totallabel],index5[i-1],0.5,0.5));
			totallabel++;
		}
	}

	  if (index!=0){
		  std::ostringstream o;
		  o << index;
		  char temp[]="xxx";
			sprintf(temp,"%d",index);

		t.push_back(new text(fd,temp,xy(0,0),1.,0.5));
	  }


	std::ofstream outfile;

	outfile.open(path);
	  pg.output(outfile);
	  outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }
  }


void print_box(boxD* thebox,const char* path,int index=0){
	double factor=2;
	double f_long=2.8;
	double f_short=1.1;
	  page pg;
	  FeynDiagram fd(pg);

//	  xy c1(-factor,factor), c2(factor,factor), c3(factor,-factor), c4(-factor,-factor);
//		 xy c[4]; c[0]=c1; c[1]=c2; c[2]=c3; c[3]=c4;
		 xy c[4]={ xy(-factor,factor), xy(factor,factor), xy(factor,-factor), xy(-factor,-factor)};

	  vertex_dot v1(fd,c[0]), v2(fd,c[1]), v3(fd,c[2]), v4(fd,c[3]);

	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][3]={"xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx"};
	  char helicities[8][2]={"1","2","3","4","5","6","7","8"};
	  xy hel_pos[8]={ xy(-f_long,-f_short),xy(-f_long,f_short), xy(-f_short,f_long), xy(f_short,f_long), xy(f_long,f_short),xy(f_long,-f_short), xy(f_short,-f_long), xy(-f_short,-f_long)};
	  for (int i=1;i<=4;i++){
		  add_line(fd,lr,thebox->l(i),c[i-1],c[(i+2)%4]);

		sprintf(helicities[2*i-2],"%s",(thebox->l(i).helicity()==-1)?"+":"-");
		sprintf(helicities[2*i-1],"%s",(thebox->l(i).helicity()==-1)?"-":"+");
		t.push_back(new text(fd,helicities[i*2-2],hel_pos[2*i-2],0.5,0.5));
		t.push_back(new text(fd,helicities[i*2-1],hel_pos[2*i-1],0.5,0.5));

	}
	  int nbr_ext;



	int totallabel=0;

	for (int l=1;l<=4;l++){
		nbr_ext=thebox->c(l).size();
		vector<xy> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double angle=5*M_PI/4.-(l*M_PI/2.)+(nbr_ext-2*i+1.)/double(nbr_ext) * M_PI/4.;
			e4[i-1]=c[l-1]+xy(2*cos(angle)*factor,2*sin(angle)*factor);
			index4[i-1]=c[l-1]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);

			add_line(fd,lr,thebox->c(l)[i-1],c[l-1],e4[i-1]);

			sprintf(temps[totallabel],"%d%s",thebox->c(l)[i-1].ind(),(thebox->c(l)[i-1].helicity()==-1)?"-":"+");
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}

	  if (index!=0){
		  std::ostringstream o;
		  o << index;
		  char temp[]="xxx";
			sprintf(temp,"%d",index);

		t.push_back(new text(fd,temp,xy(0,0),1.,0.5));
	  }


	std::ofstream outfile;

	outfile.open(path);
	  pg.output(outfile);
	  outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }
  }




void print_triangle(triangleD* thetri,const char* path,int index=0){
	double factor=2;
	double f_long=2.8;
	double f_short=1.1;
	  page pg;
	  FeynDiagram fd(pg);

	 xy c[3]={ xy(0.,factor), xy(factor,-factor), xy(-factor,-factor)};

	  vertex_dot v1(fd,c[0]), v2(fd,c[1]), v3(fd,c[2]);

	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][3]={"xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx"};
	  char helicities[6][2]={"1","2","3","4","5","6"};
	  xy hel_pos[6]={ xy(-0.9*f_long,-0.8*f_short),xy(-f_short,factor), xy(f_short,factor), xy(0.9*f_long,-0.8*f_short), xy(f_short,-f_long),xy(-f_short,-f_long)};
	  for (int i=1;i<=3;i++){

			add_line(fd,lr,thetri->l(i),c[i-1],c[(i+1)%3]);

//		  switch (thetri->l(i).t()) {
//			case g: lr.push_back(new line_spring(fd,c[(i+1)%3],c[i-1]));break;
//			case qb: lr.push_back(new line_plain(fd,c[(i+1)%3],c[i-1]));break;
//			case q:lr.push_back(new line_plain(fd,c[i-1],c[(i+1)%3]));break;
//			case qb2: lr.push_back(new line_doubleplain(fd,c[(i+1)%3],c[i-1]));break;
//			case q2:lr.push_back(new line_doubleplain(fd,c[i-1],c[(i+1)%3]));break;
//			case BH::l: lr.push_back(new line_plain(fd,c[i-1],c[(i+1)%3])); lr.back()->dashon.settrue(); lr.back()->arrowon.setfalse(); break;
//		}
		sprintf(helicities[2*i-2],"%s",(thetri->l(i).helicity()==-1)?"+":"-");
		sprintf(helicities[2*i-1],"%s",(thetri->l(i).helicity()==-1)?"-":"+");
		t.push_back(new text(fd,helicities[i*2-2],hel_pos[2*i-2],0.5,0.5));
		t.push_back(new text(fd,helicities[i*2-1],hel_pos[2*i-1],0.5,0.5));

	}
	int nbr_ext;



	int totallabel=0;

	for (int l=1;l<=3;l++){
		nbr_ext=thetri->c(l).size();
		vector<xy> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double angle=7*M_PI/6.-(l*2*M_PI/3.)+(nbr_ext-2*i+1.)/double(nbr_ext) * M_PI/3.;
			e4[i-1]=c[l-1]+xy(2*cos(angle)*factor,2*sin(angle)*factor);
			index4[i-1]=c[l-1]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);
			add_line(fd,lr,thetri->c(l)[i-1],c[l-1],e4[i-1]);
//		switch (thetri->c(l)[i-1].t()) {
//				case g: lr.push_back(new line_spring(fd,e4[i-1],c[l-1])); break;
//				case qb: lr.push_back(new line_plain(fd,e4[i-1],c[l-1])); break;
//				case q:lr.push_back(new line_plain(fd,c[l-1],e4[i-1])); break;
//				case qb2: lr.push_back(new line_doubleplain(fd,e4[i-1],c[l-1])); break;
//				case q2:lr.push_back(new line_doubleplain(fd,c[l-1],e4[i-1])); break;
//				case BH::l: lr.push_back(new line_plain(fd,c[l-1],e4[i-1])); lr.back()->dashon.settrue(); lr.back()->arrowon.setfalse(); break;
//			}
			sprintf(temps[totallabel],"%d%s",thetri->c(l)[i-1].ind(),(thetri->c(l)[i-1].helicity()==-1)?"-":"+");
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}
	  if (index!=0){
		  std::ostringstream o;
		  o << index;
		  char temp[]="xxx";
			sprintf(temp,"%d",index);
		t.push_back(new text(fd,temp,xy(0.5,-0.5),1.,0.5));
	  }


	std::ofstream outfile;

	outfile.open(path);
	  pg.output(outfile);
	  outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }
  }

void print_triangle(raw_triangle* thetri,const char* path,int index=0){
	double factor=2;
	double f_long=2.8;
	double f_short=1.1;
	  page pg;
	  FeynDiagram fd(pg);

	 xy c[3]={ xy(0.,factor), xy(factor,-factor), xy(-factor,-factor)};

	  vertex_dot v1(fd,c[0]), v2(fd,c[1]), v3(fd,c[2]);

	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][3]={"x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x"};
	  char helicities[6][2]={"1","2","3","4","5","6"};
	  xy hel_pos[6]={ xy(-0.9*f_long,-0.8*f_short),xy(-f_short,factor), xy(f_short,factor), xy(0.9*f_long,-0.8*f_short), xy(f_short,-f_long),xy(-f_short,-f_long)};
	  for (int i=1;i<=3;i++){

			add_line(fd,lr,lp,c[i-1],c[(i+1)%3]);
	}
	int nbr_ext;



	int totallabel=0;

	for (int l=1;l<=3;l++){
		nbr_ext=thetri->c(l).size();
		vector<xy> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double angle=7*M_PI/6.-(l*2*M_PI/3.)+(nbr_ext-2*i+1.)/double(nbr_ext) * M_PI/3.;
			e4[i-1]=c[l-1]+xy(2*cos(angle)*factor,2*sin(angle)*factor);
			index4[i-1]=c[l-1]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);
			add_line(fd,lr,thetri->c(l)[i-1],c[l-1],e4[i-1]);
			sprintf(temps[totallabel],"%d",thetri->indices(l)[i-1]);
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}
	  if (index!=0){
		  std::ostringstream o;
		  o << index;
		  char temp[]="xxx";
			sprintf(temp,"%d",index);
		t.push_back(new text(fd,temp,xy(0,0),1.,0.5));
	  }


	std::ofstream outfile;

	outfile.open(path);
	  pg.output(outfile);
	  outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }
  }


void print_bubble(bubbleD* thetri,const char* path,int index=0){
	double factor=2;
	double f_long=2.5;
	double f_short=2.1;
	  page pg;
	  FeynDiagram fd(pg);

	 xy c[2]={ xy(-1.5*factor,0.), xy(1.5*factor,0.)};

	  vertex_dot v1(fd,c[0]), v2(fd,c[1]);

	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][3]={"xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx"};
	  char helicities[4][2]={"1","2","3","4"};
	  xy hel_pos[4]={ xy(f_long,-f_short),xy(-f_long,-f_short), xy(-f_long,f_short), xy(f_long,f_short)};
	  for (int i=1;i<=2;i++){
			add_line(fd,lr,thetri->l(i),c[i-1],c[(i)%2]);

		lr.back()->arcthru(xy(0.,(2*(i-1.5)*factor)));
		sprintf(helicities[2*i-2],"%s",(thetri->l(i).helicity()==-1)?"+":"-");
		sprintf(helicities[2*i-1],"%s",(thetri->l(i).helicity()==-1)?"-":"+");
		t.push_back(new text(fd,helicities[i*2-2],hel_pos[2*i-2],0.5,0.5));
		t.push_back(new text(fd,helicities[i*2-1],hel_pos[2*i-1],0.5,0.5));

	}
	int nbr_ext;



	int totallabel=0;

	for (int l=1;l<=2;l++){
		nbr_ext=thetri->c(l).size();
		vector<xy> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double angle=2*M_PI-(l*M_PI)+(nbr_ext-2*i+1.)/double(nbr_ext) * 2*M_PI/3.;
			e4[i-1]=c[l-1]+xy(2*cos(angle)*factor,2*sin(angle)*factor);
			index4[i-1]=c[l-1]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);
			add_line(fd,lr,thetri->c(l)[i-1],c[l-1],e4[i-1]);
//			switch (thetri->c(l)[i-1].t()) {
//				case g: lr.push_back(new line_spring(fd,e4[i-1],c[l-1])); break;
//				case qb: lr.push_back(new line_plain(fd,e4[i-1],c[l-1])); break;
//				case q:lr.push_back(new line_plain(fd,c[l-1],e4[i-1])); break;
//				case qb2: lr.push_back(new line_doubleplain(fd,e4[i-1],c[l-1])); break;
//				case q2:lr.push_back(new line_doubleplain(fd,c[l-1],e4[i-1])); break;
//				case BH::l: lr.push_back(new line_plain(fd,c[l-1],e4[i-1])); lr.back()->dashon.settrue(); lr.back()->arrowon.setfalse(); break;
//		}
			sprintf(temps[totallabel],"%d%s",thetri->c(l)[i-1].ind(),(thetri->c(l)[i-1].helicity()==-1)?"-":"+");
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}

	  if (index!=0){
		  std::ostringstream o;
		  o << index;
		  char temp[]="xxx";
			sprintf(temp,"%d",index);
		t.push_back(new text(fd,temp,xy(0,0),1.,1.));
	  }

	std::ofstream outfile;

	outfile.open(path);
	  pg.output(outfile);
	  outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }
  }

void print_bubble(raw_bubble* thetri,const char* path,int index=0){
	double factor=2;
	double f_long=2.5;
	double f_short=2.1;
	  page pg;
	  FeynDiagram fd(pg);

	 xy c[2]={ xy(-1.5*factor,0.), xy(1.5*factor,0.)};

	  vertex_dot v1(fd,c[0]), v2(fd,c[1]);

	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][3]={"x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x"};
	  char helicities[4][2]={"1","2","3","4"};
	  xy hel_pos[4]={ xy(f_long,-f_short),xy(-f_long,-f_short), xy(-f_long,f_short), xy(f_long,f_short)};
	  for (int i=1;i<=2;i++){
			add_line(fd,lr,lp,c[i-1],c[(i)%2]);

		lr.back()->arcthru(xy(0.,(2*(i-1.5)*factor)));

	}
	int nbr_ext;



	int totallabel=0;

	for (int l=1;l<=2;l++){
		nbr_ext=thetri->c(l).size();
		vector<xy> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double angle=2*M_PI-(l*M_PI)+(nbr_ext-2*i+1.)/double(nbr_ext) * 2*M_PI/3.;
			e4[i-1]=c[l-1]+xy(2*cos(angle)*factor,2*sin(angle)*factor);
			index4[i-1]=c[l-1]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);
			add_line(fd,lr,thetri->c(l)[i-1],c[l-1],e4[i-1]);
			sprintf(temps[totallabel],"%d",thetri->indices(l)[i-1]);
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}

	  if (index!=0){
		  std::ostringstream o;
		  o << index;
		  char temp[]="xxx";
			sprintf(temp,"%d",index);
		t.push_back(new text(fd,temp,xy(0,0),1.,1.));
	  }

	std::ofstream outfile;

	outfile.open(path);
	  pg.output(outfile);
	  outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }
  }


std::string cutD_label(cutD* cd){
	char res[120];
	char props_hel[120];
	char props_type[120];

	for (int i=1;i<=cd->nc();i++){
		if ((*cd->l(i).type()) == scalar) {props_type[i-1]='a';}
		else if ((*cd->l(i).type()) == gluon) {props_type[i-1]='g';}
		else if ((*cd->l(i).type()) == photon) {props_type[i-1]='y';}
		else if ((*cd->l(i).type()) == gluino) {
			if (cd->l(i).is_anti_particle()) {
				props_type[i-1]='A';
			} else {
				props_type[i-1]='G';
			}
		}
		else if ((*cd->l(i).type()) == quark || (*cd->l(i).type()) == quark_massive) {
			if (cd->l(i).is_anti_particle()) {
				switch (cd->l(i).flavor()) {
				case 1: case -1: props_type[i-1]='q';break;
				case 2: case -2: props_type[i-1]='u';break;
				case 3: case -3: props_type[i-1]='d';break;
				case 4: case -4:props_type[i-1]='s';break;
				case 5: case -5:props_type[i-1]='c';break;
				default:props_type[i-1]='v'+cd->l(i).flavor()%5;break;
				}
			} else {
				switch (cd->l(i).flavor()) {
				case 1: case -1:props_type[i-1]='Q';break;
				case 2: case -2:props_type[i-1]='U';break;
				case 3: case -3:props_type[i-1]='D';break;
				case 4: case -4:props_type[i-1]='S';break;
				case 5: case -5:props_type[i-1]='C';break;
				default :props_type[i-1]='V'+cd->l(i).flavor()%5;break;
				}
			}
		}
		else if ((*cd->l(i).type()) == lepton) {props_type[i-1]='l';}
		else {
			_WARNING2("Untreated type in cutD_label: ",cd->l(i).type()->name());
			props_type[i-1]='b';
		}
		props_hel[i-1]=(cd->l(i).helicity()==-1)?'m':'p';
	}
	props_type[cd->nc()]='\0';
	props_hel[cd->nc()]='\0';
	sprintf(res,"%s%s%d_%ld",props_type,props_hel,cd->get_code(),cutD_external_code(*cd));
return std::string(res);

}

std::string cutD_label(raw_part* cd){
	char res[120];
	char props_hel[120];
	char props_type[120];

	for (int i=1;i<=cd->nc();i++){

	}
	sprintf(res,"g_%dU%ld",cd->get_code(),raw_part_external_code(*cd));
return std::string(res);

}

typedef cut::normal_cut_part<cut::Darren_CutD_Factory> CutType;
typedef ratext::normal_ratext RatType;


template <class CT> void print_tree_graph(CT* A,const char* path){

	bool is_rat=false;
	Cut_Part_D_Dims* CPDD = dynamic_cast<Cut_Part_D_Dims*>(A);
	if ( CPDD ){
		is_rat=true;
	}
	std::string label;
	char dot_file_path[120];
	char bubble_tex_file_path[120];
	char triangle_tex_file_path[120];
	char box_tex_file_path[120];
	char pentagon_tex_file_path[120];
	sprintf(dot_file_path,"%s/graph.dot",path);
	sprintf(bubble_tex_file_path,"%s/bubbles.tex",path);
	sprintf(triangle_tex_file_path,"%s/triangles.tex",path);
	sprintf(box_tex_file_path,"%s/boxes.tex",path);
	sprintf(pentagon_tex_file_path,"%s/pentagons.tex",path);
	_MESSAGE3("Printing graph.dot file in ",dot_file_path,".");
	std::ofstream dot_file;
	dot_file.open(dot_file_path);
	std::ofstream bubble_tex_file;
	bubble_tex_file.open(bubble_tex_file_path);
	std::ofstream triangle_tex_file;
	triangle_tex_file.open(triangle_tex_file_path);
	std::ofstream box_tex_file;
	box_tex_file.open(box_tex_file_path);
	std::ofstream pentagon_tex_file;
	pentagon_tex_file.open(pentagon_tex_file_path);

	dot_file << "digraph amplitude {" << endl ;
	std::string all_pictures;
	std::string all_html;
	std::string all_graph_pictures;
	char graph_path[120];

	dot_file << " { rank=min" << endl;

	bubble_tex_file << "\\documentclass{paper}" << endl;
	bubble_tex_file << "\\usepackage{graphicx}" << endl ;
	bubble_tex_file << "\\begin{document}" << endl ;

	triangle_tex_file << "\\documentclass{paper}" << endl;
	triangle_tex_file << "\\usepackage{graphicx}" << endl ;
	triangle_tex_file << "\\begin{document}" << endl ;

	box_tex_file << "\\documentclass{paper}" << endl;
	box_tex_file << "\\usepackage{graphicx}" << endl ;
	box_tex_file << "\\begin{document}" << endl ;

	pentagon_tex_file << "\\documentclass{paper}" << endl;
	pentagon_tex_file << "\\usepackage{graphicx}" << endl ;
	pentagon_tex_file << "\\begin{document}" << endl ;

	for (int i=0;i<A->nbr_bubbles();i++){
		label=cutD_label(A->bubble(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_bubble(A->bubble(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
		print_bubble_graph(A->bubble(i+1),path);

		bubble_tex_file << "\\includegraphics[scale=0.3]{" <<label.c_str() << ".ps}" << endl ;


	}
	dot_file << " }" << endl;
	dot_file << " { rank=same" << endl;
	for (int i=0;i<A->nbr_triangles();i++){
		label=cutD_label(A->triangle(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_triangle(A->triangle(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
		print_triangle_graph(A->triangle(i+1),path);
		triangle_tex_file << "\\includegraphics[scale=0.3]{" <<label.c_str() << ".ps}" << endl ;
	}
	dot_file << " }" << endl;
	if (is_rat){
		dot_file << " { rank=same" << endl;
	} else {
		dot_file << " { rank=max" << endl;
	}
	for (int i=0;i<A->nbr_boxes();i++){
		label=cutD_label(A->box(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_box(A->box(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
		print_box_graph(A->box(i+1),path);
		box_tex_file << "\\includegraphics[scale=0.3]{" <<label.c_str() << ".ps}" << endl ;
	}
	dot_file << " }" << endl;

	if (is_rat){
		dot_file << " { rank=max" << endl;
		for (int i=0;i<CPDD->nbr_pentagons();i++){
			label=cutD_label(CPDD->pentagon(i+1));
			dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
			sprintf(graph_path,"%s/%s.ps",path,label.c_str());
			print_pentagon(CPDD->pentagon(i+1),graph_path,i+1);
			all_pictures+=label;
			all_pictures+=".png ";
			all_pictures+=label;
			all_pictures+="_small.png ";
			all_html+=label;
			all_html+=".html ";
			all_graph_pictures+=label;
			all_graph_pictures+="_graph.png ";
//			print_box_graph(CPDD->pentagon(i+1),path);
			pentagon_tex_file << "\\includegraphics[scale=0.3]{" <<label.c_str() << ".ps}" << endl ;
		}
		dot_file << " }" << endl;
	}

	for (int i=0;i<A->nbr_bubbles();i++){
		std::string bub_label=cutD_label(A->bubble(i+1));
		for (int j=0;j<A->bubble(i+1)->daughters_nbr();j++){
			dot_file << bub_label << " -> " << cutD_label(A->bubble(i+1)->get_daughter(j+1)) << ";" <<endl;
		}
	}


	for (int i=0;i<A->nbr_triangles();i++){
		std::string tri_label=cutD_label(A->triangle(i+1));
		for (int j=0;j<A->triangle(i+1)->daughters_nbr();j++){
			dot_file << cutD_label(A->triangle(i+1)) << " -> " << cutD_label(A->triangle(i+1)->get_daughter(j+1)) << ";" <<endl;
		}
	}

	for (int i=0;i<A->nbr_triangles();i++){
		std::string tri_label=cutD_label(A->triangle(i+1));
		for (int j=0;j<A->triangle(i+1)->parents_nbr();j++){
			dot_file << cutD_label(A->triangle(i+1)) << " -> " << cutD_label(A->triangle(i+1)->get_parent(j+1)) << ";" <<endl;
		}
	}

	for (int i=0;i<A->nbr_boxes();i++){
		std::string bub_label=cutD_label(A->box(i+1));
		for (int j=0;j<A->box(i+1)->parents_nbr();j++){
			dot_file << bub_label << " -> " << cutD_label(A->box(i+1)->get_parent(j+1)) << ";" <<endl;
		}
	}

	if (is_rat){

		for (int i=0;i<CPDD->nbr_boxes();i++){
			std::string tri_label=cutD_label(CPDD->box(i+1));
			boxD_D_Dims* bd=dynamic_cast<boxD_D_Dims*>(CPDD->box(i+1));
			for (int j=0;j<bd->daughters_nbr();j++){
				dot_file << cutD_label(CPDD->box(i+1)) << " -> " << cutD_label(bd->get_daughter(j+1)) << ";" <<endl;
			}
		}

		for (int i=0;i<CPDD->nbr_pentagons();i++){
			std::string bub_label=cutD_label(CPDD->pentagon(i+1));
			for (int j=0;j<CPDD->pentagon(i+1)->parents_nbr();j++){
				dot_file << bub_label << " -> " << cutD_label(CPDD->pentagon(i+1)->get_parent(j+1)) << ";" <<endl;
			}
		}

	}

dot_file << "}"<<endl;
dot_file.close();
bubble_tex_file << "\\end{document}" << endl ;
bubble_tex_file.close();
triangle_tex_file << "\\end{document}" << endl ;
triangle_tex_file.close();
box_tex_file << "\\end{document}" << endl ;
box_tex_file.close();

	std::ofstream makefile;
	char makefile_path[120];
	sprintf(makefile_path,"%s/makefile",path);
	_MESSAGE3("Printing makefile in ",makefile_path,".");
	makefile.open(makefile_path);

makefile << "%.epsi:\t%.ps" << endl;
makefile << "\tps2epsi $<" << endl;
makefile << "%.png:\t%.epsi" << endl;
makefile << "\tconvert $< $@" << endl;
makefile << "%_small.png:\t%.png" << endl;
makefile << "\tconvert  -scale 200x200  $< $@" << endl;
makefile << "graph.png:\tgraph.dot " << all_pictures  << endl;
makefile << "\tdot -Tpng -o graph.png graph.dot  "  << endl;
makefile << "graph.imap:\tgraph.dot "   << all_pictures<< endl;
makefile << "\tdot -Timap -o graph.imap graph.dot  "  << endl;
makefile << "graph.cmapx:\tgraph.dot"  << endl;
makefile << "\tdot -Tcmap -o graph.cmapx graph.dot  "  << endl;

makefile << "%_graph.png:\t%.dot " << all_pictures  << endl;
makefile << "\tdot -Tpng -o $*_graph.png $*.dot  "  << endl;
makefile << "%.imap:\t%.dot "   << all_pictures<< endl;
makefile << "\tdot -Timap -o $*.imap $*.dot  "  << endl;
makefile << "%.cmapx:\t%.dot"  << endl;
makefile << "\tdot -Tcmap -o $*.cmapx $*.dot  "  << endl;


makefile << "index.html:\tgraph.png graph.imap graph.cmapx"  << endl;
makefile << "\tcat index.html.in1 > index.html  "  << endl;

makefile << "\tebb graph.png"  << endl;
makefile << "\techo -n \"width: \" >>  index.html" << endl;
makefile << "\techo -n `cat graph.bb | grep BoundingBox | grep -m 1 -o \"[0-9]*[0-9][0-9] [1-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>index.html" << endl;
makefile << "\techo -n \"; height: \" >>  index.html " << endl;
makefile << "\techo -n `cat graph.bb | grep BoundingBox | grep -m 1 -o \"[1-9] [0-9]*[0-9][0-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>index.html " << endl;
makefile << "\techo  \";\" >>  index.html" << endl;
makefile << "\techo \"  graph.png\"  >>  index.html " << endl;

makefile << "\tcat index.html.in2 >> index.html  "  << endl;
makefile << "\tcat graph.cmapx >> index.html  "  << endl;
makefile << "\tcat index.html.in3 >> index.html  "  << endl;

makefile << "%.html:\t%.imap %.cmapx %_graph.png"  << endl;
makefile << "\tcat $*.html.in1 > $@  "  << endl;

makefile << "\tebb $*_graph.png"  << endl;
makefile << "\techo -n \"<img style=\\\"width: \">>  $@ "  << endl;
makefile << "\techo -n `cat $*_graph.bb | grep BoundingBox | grep -m 1 -o \"[0-9]*[0-9][0-9] [1-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>$@" << endl;
makefile << "\techo -n \"; height: \" >>  $@ " << endl;
makefile << "\techo -n `cat $*_graph.bb | grep BoundingBox | grep -m 1 -o \"[1-9] [0-9]*[0-9][0-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>$@ " << endl;
makefile << "\techo  \";\" >>  $@ " << endl;
makefile << "\techo \"  $*_graph.png\"  >>  $@ " << endl;


//makefile << "\tpnginfo -f \"<img style=\\\"width: \\wpx; height: \\hpx;\\\"\\n\"  $*_graph.png >>  $@  "  << endl;
makefile << "\tcat $*.html.in2 >> $@  "  << endl;
makefile << "\tcat $*.cmapx >> $@  "  << endl;
makefile << "\tcat $*.html.in3 >> $@  "  << endl;


makefile << "all:\tindex.html " << all_pictures << all_graph_pictures << all_html << endl;
makefile.close();

	std::ofstream html;
	char html_path[120];
	sprintf(html_path,"%s/index.html.in1",path);
	_MESSAGE3("Generating index.html.in in ",makefile_path,".");
	html.open(html_path);

	html<< "<html>" <<endl;
	html<< "<head>" <<endl;
	html<< "  <meta content=\"text/html; charset=ISO-8859-1\"" <<endl;
	html<< " http-equiv=\"content-type\">" <<endl;
	html<< "  <title>graph</title>" <<endl;
	html<< "</head>" <<endl;
	html<< "<body>" <<endl;
	html<< "<img style=\""<<endl;
	html.close();
	sprintf(html_path,"%s/index.html.in2",path);
	html.open(html_path);
	//	width: 5607px; height: 935px; provided by the makefile script
	html<< "\" usemap=\"#amplitude\"" <<endl;
	html<< " alt=\"graph\"" <<endl;
	html<< " src=\"graph.png\">" <<endl;
	html<< "<map id=\"amplitude\" name=\"amplitude\">" <<endl;

	html.close();
	sprintf(html_path,"%s/index.html.in3",path);
	html.open(html_path);
	html<< "</body>" <<endl;
	html<< "</html>" <<endl;
	html.close();
}

void print_tree_graph_non_zero(Cut_Part* A,const char* path,mom_conf& mc,const vector<int>& ind){
	std::string label;
	char dot_file_path[120];
	char bubble_tex_file_path[120];
	char triangle_tex_file_path[120];
	char box_tex_file_path[120];
	sprintf(dot_file_path,"%s/graph.dot",path);
	sprintf(bubble_tex_file_path,"%s/bubbles.tex",path);
	sprintf(triangle_tex_file_path,"%s/triangles.tex",path);
	sprintf(box_tex_file_path,"%s/boxes.tex",path);
	_MESSAGE3("Printing graph.dot file in ",dot_file_path,".");
	std::ofstream dot_file;
	dot_file.open(dot_file_path);
	std::ofstream bubble_tex_file;
	bubble_tex_file.open(bubble_tex_file_path);
	std::ofstream triangle_tex_file;
	triangle_tex_file.open(triangle_tex_file_path);
	std::ofstream box_tex_file;
	box_tex_file.open(box_tex_file_path);

	dot_file << "digraph amplitude {" << endl ;
	std::string all_pictures;
	std::string all_html;
	std::string all_graph_pictures;
	char graph_path[120];

	dot_file << " { rank=min" << endl;

	bubble_tex_file << "\\documentclass{paper}" << endl;
	bubble_tex_file << "\\usepackage{graphicx}" << endl ;
	bubble_tex_file << "\\begin{document}" << endl ;

	triangle_tex_file << "\\documentclass{paper}" << endl;
	triangle_tex_file << "\\usepackage{graphicx}" << endl ;
	triangle_tex_file << "\\begin{document}" << endl ;

	box_tex_file << "\\documentclass{paper}" << endl;
	box_tex_file << "\\usepackage{graphicx}" << endl ;
	box_tex_file << "\\begin{document}" << endl ;

	for (int i=0;i<A->nbr_bubbles();i++){
		if ( abs(A->bubble(i+1)->get_value(mc,ind)) > 1e-7 ){
		label=cutD_label(A->bubble(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_bubble(A->bubble(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
		print_bubble_graph(A->bubble(i+1),path);

		bubble_tex_file << "\\includegraphics[scale=0.3]{" <<label.c_str() << ".ps}" << endl ;
		}

	}
	dot_file << " }" << endl;
	dot_file << " { rank=same" << endl;
	for (int i=0;i<A->nbr_triangles();i++){
		if (abs (A->triangle(i+1)->get_value(mc,ind))> 1e-7){
		label=cutD_label(A->triangle(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_triangle(A->triangle(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
		print_triangle_graph(A->triangle(i+1),path);
		triangle_tex_file << "\\includegraphics[scale=0.3]{" <<label.c_str() << ".ps}" << endl ;
		}
		}
	dot_file << " }" << endl;
	dot_file << " { rank=max" << endl;
	for (int i=0;i<A->nbr_boxes();i++){
		if (abs(A->box(i+1)->get_value(mc,ind))>1e-7){
		label=cutD_label(A->box(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_box(A->box(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
		print_box_graph(A->box(i+1),path);
		box_tex_file << "\\includegraphics[scale=0.3]{" <<label.c_str() << ".ps}" << endl ;
		}
	}
	dot_file << " }" << endl;

	for (int i=0;i<A->nbr_bubbles();i++){
		std::string bub_label=cutD_label(A->bubble(i+1));
		if (abs (A->bubble(i+1)->get_value(mc,ind))> 1e-7){
			for (int j=0;j<A->bubble(i+1)->daughters_nbr();j++){
				if (abs (A->bubble(i+1)->get_daughter(j+1)->get_value(mc,ind))> 1e-7){
						dot_file << bub_label << " -> " << cutD_label(A->bubble(i+1)->get_daughter(j+1)) << ";" <<endl;
				}
			}
		}
	}


	for (int i=0;i<A->nbr_triangles();i++){
		std::string tri_label=cutD_label(A->triangle(i+1));
		if (abs (A->triangle(i+1)->get_value(mc,ind))> 1e-7){
			for (int j=0;j<A->triangle(i+1)->daughters_nbr();j++){
				if (abs (A->triangle(i+1)->get_daughter(j+1)->get_value(mc,ind))> 1e-7){
					dot_file << cutD_label(A->triangle(i+1)) << " -> " << cutD_label(A->triangle(i+1)->get_daughter(j+1)) << ";" <<endl;
				}
			}
		}
	}

	for (int i=0;i<A->nbr_triangles();i++){
		if (abs (A->triangle(i+1)->get_value(mc,ind))> 1e-7){
			std::string tri_label=cutD_label(A->triangle(i+1));
			for (int j=0;j<A->triangle(i+1)->parents_nbr();j++){
				if (abs (A->triangle(i+1)->get_parent(j+1)->get_value(mc,ind))> 1e-7){
					dot_file << cutD_label(A->triangle(i+1)) << " -> " << cutD_label(A->triangle(i+1)->get_parent(j+1)) << ";" <<endl;
				}
			}
		}
	}

	for (int i=0;i<A->nbr_boxes();i++){
		if (abs (A->box(i+1)->get_value(mc,ind))> 1e-7){
			std::string bub_label=cutD_label(A->box(i+1));
			for (int j=0;j<A->box(i+1)->parents_nbr();j++){
				if (abs (A->box(i+1)->get_parent(j+1)->get_value(mc,ind))> 1e-7){
					dot_file << bub_label << " -> " << cutD_label(A->box(i+1)->get_parent(j+1)) << ";" <<endl;
				}
			}
		}
	}

dot_file << "}"<<endl;
dot_file.close();
bubble_tex_file << "\\end{document}" << endl ;
bubble_tex_file.close();
triangle_tex_file << "\\end{document}" << endl ;
triangle_tex_file.close();
box_tex_file << "\\end{document}" << endl ;
box_tex_file.close();

	std::ofstream makefile;
	char makefile_path[120];
	sprintf(makefile_path,"%s/makefile",path);
	_MESSAGE3("Printing makefile in ",makefile_path,".");
	makefile.open(makefile_path);

makefile << "%.epsi:\t%.ps" << endl;
makefile << "\tps2epsi $<" << endl;
makefile << "%.png:\t%.epsi" << endl;
makefile << "\tconvert $< $@" << endl;
makefile << "%_small.png:\t%.png" << endl;
makefile << "\tconvert  -scale 200x200  $< $@" << endl;
makefile << "graph.png:\tgraph.dot " << all_pictures  << endl;
makefile << "\tdot -Tpng -o graph.png graph.dot  "  << endl;
makefile << "graph.imap:\tgraph.dot "   << all_pictures<< endl;
makefile << "\tdot -Timap -o graph.imap graph.dot  "  << endl;
makefile << "graph.cmapx:\tgraph.dot"  << endl;
makefile << "\tdot -Tcmap -o graph.cmapx graph.dot  "  << endl;

makefile << "%_graph.png:\t%.dot " << all_pictures  << endl;
makefile << "\tdot -Tpng -o $*_graph.png $*.dot  "  << endl;
makefile << "%.imap:\t%.dot "   << all_pictures<< endl;
makefile << "\tdot -Timap -o $*.imap $*.dot  "  << endl;
makefile << "%.cmapx:\t%.dot"  << endl;
makefile << "\tdot -Tcmap -o $*.cmapx $*.dot  "  << endl;


makefile << "index.html:\tgraph.png graph.imap graph.cmapx"  << endl;
makefile << "\tcat index.html.in1 > index.html  "  << endl;

makefile << "\tebb graph.png"  << endl;
makefile << "\techo -n \"width: \" >>  index.html" << endl;
makefile << "\techo -n `cat graph.bb | grep BoundingBox | grep -m 1 -o \"[0-9]*[0-9][0-9] [1-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>index.html" << endl;
makefile << "\techo -n \"; height: \" >>  index.html " << endl;
makefile << "\techo -n `cat graph.bb | grep BoundingBox | grep -m 1 -o \"[1-9] [0-9]*[0-9][0-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>index.html " << endl;
makefile << "\techo  \";\" >>  index.html" << endl;
makefile << "\techo \"  graph.png\"  >>  index.html " << endl;

makefile << "\tcat index.html.in2 >> index.html  "  << endl;
makefile << "\tcat graph.cmapx >> index.html  "  << endl;
makefile << "\tcat index.html.in3 >> index.html  "  << endl;

makefile << "%.html:\t%.imap %.cmapx %_graph.png"  << endl;
makefile << "\tcat $*.html.in1 > $@  "  << endl;

makefile << "\tebb $*_graph.png"  << endl;
makefile << "\techo -n \"<img style=\\\"width: \">>  $@ "  << endl;
makefile << "\techo -n `cat $*_graph.bb | grep BoundingBox | grep -m 1 -o \"[0-9]*[0-9][0-9] [1-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>$@" << endl;
makefile << "\techo -n \"; height: \" >>  $@ " << endl;
makefile << "\techo -n `cat $*_graph.bb | grep BoundingBox | grep -m 1 -o \"[1-9] [0-9]*[0-9][0-9]\" | grep -o \"[0-9]*[0-9][0-9]\" ` >>$@ " << endl;
makefile << "\techo  \";\" >>  $@ " << endl;
makefile << "\techo \"  $*_graph.png\"  >>  $@ " << endl;


//makefile << "\tpnginfo -f \"<img style=\\\"width: \\wpx; height: \\hpx;\\\"\\n\"  $*_graph.png >>  $@  "  << endl;
makefile << "\tcat $*.html.in2 >> $@  "  << endl;
makefile << "\tcat $*.cmapx >> $@  "  << endl;
makefile << "\tcat $*.html.in3 >> $@  "  << endl;


makefile << "all:\tindex.html " << all_pictures << all_graph_pictures << all_html << endl;
makefile.close();

	std::ofstream html;
	char html_path[120];
	sprintf(html_path,"%s/index.html.in1",path);
	_MESSAGE3("Generating index.html.in in ",makefile_path,".");
	html.open(html_path);

	html<< "<html>" <<endl;
	html<< "<head>" <<endl;
	html<< "  <meta content=\"text/html; charset=ISO-8859-1\"" <<endl;
	html<< " http-equiv=\"content-type\">" <<endl;
	html<< "  <title>graph</title>" <<endl;
	html<< "</head>" <<endl;
	html<< "<body>" <<endl;
	html<< "<img style=\""<<endl;
	html.close();
	sprintf(html_path,"%s/index.html.in2",path);
	html.open(html_path);
	//	width: 5607px; height: 935px; provided by the makefile script
	html<< "\" usemap=\"#amplitude\"" <<endl;
	html<< " alt=\"graph\"" <<endl;
	html<< " src=\"graph.png\">" <<endl;
	html<< "<map id=\"amplitude\" name=\"amplitude\">" <<endl;

	html.close();
	sprintf(html_path,"%s/index.html.in3",path);
	html.open(html_path);
	html<< "</body>" <<endl;
	html<< "</html>" <<endl;
	html.close();
}

void print_triangle_graph(triangleD* cd,const char* path){
	std::string label;
	std::string this_label=cutD_label(cd);

	char dot_file_path[120];
	sprintf(dot_file_path,"%s/%s.dot",path,this_label.c_str());
	std::ofstream dot_file;
	dot_file.open(dot_file_path);

	dot_file << "digraph " << this_label << " {" << endl ;

	std::string all_pictures;
	char graph_path[120];
	for (int i=0;i<cd->parents_nbr();i++){
		label=cutD_label(cd->get_parent(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<".png\"];"<<std::endl;
	}

	dot_file << this_label << " [label=\"\",href=\"" <<this_label<<".png\",shapefile=\"" << this_label <<".png\"];"<<std::endl;


	for (int i=0;i<cd->daughters_nbr();i++){
		label=cutD_label(cd->get_daughter(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<".png\"];"<<std::endl;
	}


	for (int j=0;j<cd->daughters_nbr();j++){
		dot_file << this_label << " -> " << cutD_label(cd->get_daughter(j+1)) << ";" <<endl;
	}

	for (int j=0;j<cd->parents_nbr();j++){
		dot_file << cutD_label(cd->get_parent(j+1)) << " -> " << this_label << ";" <<endl;
	}

	dot_file << "}"<<endl;
	dot_file.close();

	std::ofstream html;
	char html_path[120];
	sprintf(html_path,"%s/%s.html.in1",path,this_label.c_str());
	html.open(html_path);

	html<< "<html>" <<endl;
	html<< "<head>" <<endl;
	html<< "  <meta content=\"text/html; charset=ISO-8859-1\"" <<endl;
	html<< " http-equiv=\"content-type\">" <<endl;
	html<< "  <title>graph</title>" <<endl;
	html<< "</head>" <<endl;
	html<< "<body>" <<endl;
//	html<< "<img style=\"width: 5607px; height: 935px;\"
	html.close();
	sprintf(html_path,"%s/%s.html.in2",path,this_label.c_str());
	html.open(html_path);
	html << "usemap=\"#" << this_label << "\"" <<endl;
	html<< " alt=\"graph\"" <<endl;
	html<< " src=\""<< this_label << "_graph.png\">" <<endl;
	html<< "<map id=\"" << this_label << "\" name=\""<< this_label << "\">" <<endl;

	html.close();

	sprintf(html_path,"%s/%s.html.in3",path,this_label.c_str());
	html.open(html_path);
	html<< "</body>" <<endl;
	html<< "</html>" <<endl;
	html.close();

}

void print_box_graph(boxD* cd,const char* path){
	std::string label;
	std::string this_label=cutD_label(cd);

	char dot_file_path[120];
	sprintf(dot_file_path,"%s/%s.dot",path,this_label.c_str());
	std::ofstream dot_file;
	dot_file.open(dot_file_path);

	dot_file << "digraph " << this_label << " {" << endl ;

	std::string all_pictures;
	char graph_path[120];
	for (int i=0;i<cd->parents_nbr();i++){
		label=cutD_label(cd->get_parent(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<".png\"];"<<std::endl;
	}

	dot_file << this_label << " [label=\"\",href=\"" <<this_label<<".html\",shapefile=\"" << this_label <<".png\"];"<<std::endl;


	for (int j=0;j<cd->parents_nbr();j++){
		dot_file << cutD_label(cd->get_parent(j+1)) << " -> " << this_label << ";" <<endl;
	}

	dot_file << "}"<<endl;
	dot_file.close();

	std::ofstream html;
	char html_path[120];
	sprintf(html_path,"%s/%s.html.in1",path,this_label.c_str());
	html.open(html_path);

	html<< "<html>" <<endl;
	html<< "<head>" <<endl;
	html<< "  <meta content=\"text/html; charset=ISO-8859-1\"" <<endl;
	html<< " http-equiv=\"content-type\">" <<endl;
	html<< "  <title>graph</title>" <<endl;
	html<< "</head>" <<endl;
	html<< "<body>" <<endl;
//	html<< "<img style=\"width: 5607px; height: 935px;\"
	html.close();
	sprintf(html_path,"%s/%s.html.in2",path,this_label.c_str());
	html.open(html_path);
	html << "usemap=\"#" << this_label << "\"" <<endl;
	html<< " alt=\"graph\"" <<endl;
	html<< " src=\""<< this_label << "_graph.png\">" <<endl;
	html<< "<map id=\"" << this_label << "\" name=\""<< this_label << "\">" <<endl;

	html.close();

	sprintf(html_path,"%s/%s.html.in3",path,this_label.c_str());
	html.open(html_path);
	html<< "</body>" <<endl;
	html<< "</html>" <<endl;
	html.close();

}



void print_bubble_graph(bubbleD* cd,const char* path){
	std::string label;
	std::string this_label=cutD_label(cd);

	char dot_file_path[120];
	sprintf(dot_file_path,"%s/%s.dot",path,this_label.c_str());
	std::ofstream dot_file;
	dot_file.open(dot_file_path);

	dot_file << "digraph " << this_label << " {" << endl ;

	std::string all_pictures;
	char graph_path[120];
	dot_file << this_label << " [label=\"\",href=\"" <<this_label<<".html\",shapefile=\"" << this_label <<".png\"];"<<std::endl;


	for (int i=0;i<cd->daughters_nbr();i++){
		label=cutD_label(cd->get_daughter(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<".png\"];"<<std::endl;
	}


	for (int j=0;j<cd->daughters_nbr();j++){
		dot_file << this_label << " -> " << cutD_label(cd->get_daughter(j+1)) << ";" <<endl;
	}


	dot_file << "}"<<endl;
	dot_file.close();

	std::ofstream html;
	char html_path[120];
	sprintf(html_path,"%s/%s.html.in1",path,this_label.c_str());
	html.open(html_path);

	html<< "<html>" <<endl;
	html<< "<head>" <<endl;
	html<< "  <meta content=\"text/html; charset=ISO-8859-1\"" <<endl;
	html<< " http-equiv=\"content-type\">" <<endl;
	html<< "  <title>graph</title>" <<endl;
	html<< "</head>" <<endl;
	html<< "<body>" <<endl;
//	html<< "<img style=\"width: 5607px; height: 935px;\"
	html.close();
	sprintf(html_path,"%s/%s.html.in2",path,this_label.c_str());
	html.open(html_path);
	html << "usemap=\"#" << this_label << "\"" <<endl;
	html<< " alt=\"graph\"" <<endl;
	html<< " src=\""<< this_label << "_graph.png\">" <<endl;
	html<< "<map id=\"" << this_label << "\" name=\""<< this_label << "\">" <<endl;

	html.close();

	sprintf(html_path,"%s/%s.html.in3",path,this_label.c_str());
	html.open(html_path);
	html<< "</body>" <<endl;
	html<< "</html>" <<endl;
	html.close();

}


void print_tree_graph(OneLoopRawAmplitude* A,const char* path){
	std::string label;
	char dot_file_path[120];
	sprintf(dot_file_path,"%s/graph.dot",path);
	_MESSAGE3("Printing graph.dot file in ",dot_file_path,".");
	std::ofstream dot_file;
	dot_file.open(dot_file_path);

	dot_file << "digraph amplitude {" << endl ;
	std::string all_pictures;
	std::string all_html;
	std::string all_graph_pictures;
	char graph_path[120];

	dot_file << " { rank=min" << endl;

	for (int i=0;i<A->nbr_bubbles();i++){
		label=cutD_label(A->bubble(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_bubble(A->bubble(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
//		print_bubble_graph(A->bubble(i+1),path);
	}
	dot_file << " }" << endl;

	dot_file << " { rank=same" << endl;
	for (int i=0;i<A->nbr_triangles();i++){
		label=cutD_label(A->triangle(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_triangle(A->triangle(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
//		print_triangle_graph(A->triangle(i+1),path);
	}
	dot_file << " }" << endl;
	dot_file << " { rank=max" << endl;
	for (int i=0;i<A->nbr_boxes();i++){
		label=cutD_label(A->box(i+1));
		dot_file << label << " [label=\"\",href=\"" <<label<<".html\",shapefile=\"" << label <<"_small.png\"];"<<std::endl;
		sprintf(graph_path,"%s/%s.ps",path,label.c_str());
		print_box(A->box(i+1),graph_path,i+1);
		all_pictures+=label;
		all_pictures+=".png ";
		all_pictures+=label;
		all_pictures+="_small.png ";
		all_html+=label;
		all_html+=".html ";
		all_graph_pictures+=label;
		all_graph_pictures+="_graph.png ";
//		print_box_graph(A->box(i+1),path);
	}
	dot_file << " }" << endl;

	for (int i=0;i<A->nbr_bubbles();i++){
		std::string bub_label=cutD_label(A->bubble(i+1));
		for (int j=0;j<A->bubble(i+1)->daughters_nbr();j++){
			dot_file << bub_label << " -> " << cutD_label(A->bubble(i+1)->get_daughter(j+1)) << ";" <<endl;
		}
	}


	for (int i=0;i<A->nbr_triangles();i++){
		std::string tri_label=cutD_label(A->triangle(i+1));
		for (int j=0;j<A->triangle(i+1)->daughters_nbr();j++){
			dot_file << cutD_label(A->triangle(i+1)) << " -> " << cutD_label(A->triangle(i+1)->get_daughter(j+1)) << ";" <<endl;
		}
	}

	for (int i=0;i<A->nbr_triangles();i++){
		std::string tri_label=cutD_label(A->triangle(i+1));
		for (int j=0;j<A->triangle(i+1)->parents_nbr();j++){
			dot_file << cutD_label(A->triangle(i+1)) << " -> " << cutD_label(A->triangle(i+1)->get_parent(j+1)) << ";" <<endl;
		}
	}

	for (int i=0;i<A->nbr_boxes();i++){
		std::string bub_label=cutD_label(A->box(i+1));
		for (int j=0;j<A->box(i+1)->parents_nbr();j++){
			dot_file << bub_label << " -> " << cutD_label(A->box(i+1)->get_parent(j+1)) << ";" <<endl;
		}
	}

dot_file << "}"<<endl;
dot_file.close();

	std::ofstream makefile;
	char makefile_path[120];
	sprintf(makefile_path,"%s/makefile",path);
	_MESSAGE3("Printing makefile in ",makefile_path,".");
	makefile.open(makefile_path);

makefile << "%.epsi:\t%.ps" << endl;
makefile << "\tps2epsi $<" << endl;
makefile << "%.png:\t%.epsi" << endl;
makefile << "\tconvert $< $@" << endl;
makefile << "%_small.png:\t%.png" << endl;
makefile << "\tconvert  -scale 200x200  $< $@" << endl;
makefile << "graph.png:\tgraph.dot " << all_pictures  << endl;
makefile << "\tdot -Tpng -o graph.png graph.dot  "  << endl;
makefile << "graph.imap:\tgraph.dot "   << all_pictures<< endl;
makefile << "\tdot -Timap -o graph.imap graph.dot  "  << endl;
makefile << "graph.cmapx:\tgraph.dot"  << endl;
makefile << "\tdot -Tcmap -o graph.cmapx graph.dot  "  << endl;

makefile << "%_graph.png:\t%.dot " << all_pictures  << endl;
makefile << "\tdot -Tpng -o $*_graph.png $*.dot  "  << endl;
makefile << "%.imap:\t%.dot "   << all_pictures<< endl;
makefile << "\tdot -Timap -o $*.imap $*.dot  "  << endl;
makefile << "%.cmapx:\t%.dot"  << endl;
makefile << "\tdot -Tcmap -o $*.cmapx $*.dot  "  << endl;


makefile << "index.html:\tgraph.png graph.imap graph.cmapx"  << endl;
makefile << "\tcat index.html.in1 > index.html  "  << endl;
makefile << "\tpnginfo -f \"width: \\wpx; height: \\hpx;\\n\"  graph.png >>  index.html  "  << endl;
makefile << "\tcat index.html.in2 >> index.html  "  << endl;
makefile << "\tcat graph.cmapx >> index.html  "  << endl;
makefile << "\tcat index.html.in3 >> index.html  "  << endl;

makefile << "%.html:\t%.imap %.cmapx %_graph.png"  << endl;
makefile << "\tcat $*.html.in1 > $@  "  << endl;
makefile << "\tpnginfo -f \"<img style=\\\"width: \\wpx; height: \\hpx;\\\"\\n\"  $*_graph.png >>  $@  "  << endl;
makefile << "\tcat $*.html.in2 >> $@  "  << endl;
makefile << "\tcat $*.cmapx >> $@  "  << endl;
makefile << "\tcat $*.html.in3 >> $@  "  << endl;


makefile << "all:\tindex.html " << all_pictures << all_graph_pictures << all_html << endl;
makefile.close();

	std::ofstream html;
	char html_path[120];
	sprintf(html_path,"%s/index.html.in1",path);
	_MESSAGE3("Generating index.html.in in ",makefile_path,".");
	html.open(html_path);

	html<< "<html>" <<endl;
	html<< "<head>" <<endl;
	html<< "  <meta content=\"text/html; charset=ISO-8859-1\"" <<endl;
	html<< " http-equiv=\"content-type\">" <<endl;
	html<< "  <title>graph</title>" <<endl;
	html<< "</head>" <<endl;
	html<< "<body>" <<endl;
	html<< "<img style=\""<<endl;
	html.close();
	sprintf(html_path,"%s/index.html.in2",path);
	html.open(html_path);
	//	width: 5607px; height: 935px; provided by the makefile script
	html<< "\" usemap=\"#amplitude\"" <<endl;
	html<< " alt=\"graph\"" <<endl;
	html<< " src=\"graph.png\">" <<endl;
	html<< "<map id=\"amplitude\" name=\"amplitude\">" <<endl;

	html.close();
	sprintf(html_path,"%s/index.html.in3",path);
	html.open(html_path);
	html<< "</body>" <<endl;
	html<< "</html>" <<endl;
	html.close();
}

void print_cut_part_graph(Cut_Part_base* CPB,const char* path){
	CutType* CP=dynamic_cast<CutType*>(CPB);
	if ( CP != 0){
		print_tree_graph(CP,path);
		return;
	}
	Cut_Part_D_Dims* RP=dynamic_cast<Cut_Part_D_Dims*>(CPB);
	if ( RP != 0){
		print_tree_graph(RP,path);
		return;
	}


	_MESSAGE3("Sorry, can't print this kind of cut part (" , typeid(*CP).name() ,").");
}



void print_cutD(const cutD& cd,const char* filename, const std::vector<int>& ind){

	double factor=2;
	double f_long=2.8;
	double f_short=1.1;
	  page pg;

	  FeynDiagram fd(pg);

	  int n=cd.nc();
//	  vertex_dot(fd,xy(1,1));


	  vector<xy*> corners;
	  vector<vertex_dot*> v;
	  double total=6.2830;

	  for (int i=0;i<n;i++){
		  double angle=-i*total/double(n);
		  corners.push_back(new xy(factor*cos(angle),factor*sin(angle)));
		  v.push_back(new vertex_dot(fd,*corners[i]));
	  };


	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][3]={"xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx"};

	  for (int i=1;i<=n;i++){
		  add_line(fd,lr,cd.l(i),*corners[i-1],*corners[(i-2+n)%n]);
		  if ( n == 2){
			  lr.back()->arcthru(xy(0.,-(2*(i-1.5)*factor)));
		  }
	  }

	  int nbr_ext;
//
//
//
	int totallabel=0;
//
	for (int l=0;l<n;l++){
		nbr_ext=cd.c(l+1).size();
		vector<xy*> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double base_angle=-(l*2.*M_PI/double(n));
			double angle=base_angle+(nbr_ext-2*i+1.)/double(nbr_ext) * M_PI/4.;
			e4[i-1]=new xy(*corners[l]+xy(2*cos(angle)*factor,2*sin(angle)*factor));
			index4[i-1]=*corners[l]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);

			add_line(fd,lr,cd.c(l+1)[i-1],*corners[l],*e4[i-1]);

			sprintf(temps[totallabel],"%d",ind[cd.c(l+1)[i-1].ind()-1]);
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}

	std::ofstream outfile;

	outfile.open(filename);
	pg.output(outfile);
	outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }

}



void print_cutD(const cutD& cd,const char* filename){

	double factor=2;
	double f_long=2.8;
	double f_short=1.1;
	  page pg;

	  FeynDiagram fd(pg);

	  int n=cd.nc();
//	  vertex_dot(fd,xy(1,1));


	  vector<xy*> corners;
	  vector<vertex_dot*> v;
	  double total=6.2830;

	  for (int i=0;i<n;i++){
		  double angle=-i*total/double(n);
		  corners.push_back(new xy(factor*cos(angle),factor*sin(angle)));
		  v.push_back(new vertex_dot(fd,*corners[i]));
	  };


	  vector<line_raw*> lr(0);
	  vector<text*> t;
	  char temps[20][3]={"xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx","xx"};

	  for (int i=1;i<=n;i++){
		  add_line(fd,lr,cd.l(i),*corners[i-1],*corners[(i-2+n)%n]);
		  if ( n == 2){
			  lr.back()->arcthru(xy(0.,-(2*(i-1.5)*factor)));
		  }
	  }

	  int nbr_ext;
//
//
//
	int totallabel=0;
//
	for (int l=0;l<n;l++){
		nbr_ext=cd.c(l+1).size();
		vector<xy*> e4(nbr_ext);
		vector<xy> index4(nbr_ext);
		for (int i=1;i<=nbr_ext;i++){
			double base_angle=-(l*2.*M_PI/double(n));
			double angle=base_angle+(nbr_ext-2*i+1.)/double(nbr_ext) * M_PI/4.;
			e4[i-1]=new xy(*corners[l]+xy(2*cos(angle)*factor,2*sin(angle)*factor));
			index4[i-1]=*corners[l]+xy(2.5*cos(angle)*factor,2.5*sin(angle)*factor);

			add_line(fd,lr,cd.c(l+1)[i-1],*corners[l],*e4[i-1]);

			sprintf(temps[totallabel],"%d",cd.c(l+1)[i-1].ind());
			t.push_back(new text(fd,temps[totallabel],index4[i-1],0.5,0.5));
			totallabel++;
		}
	}

	std::ofstream outfile;

	outfile.open(filename);
	pg.output(outfile);
	outfile.close();
	  for (int k=0;k<lr.size();k++) {
		  delete lr[k];
	  }
	  for (int k=0;k<t.size();k++) {
		  delete t[k];
	  }

}

void display_cut(const cutD& cd){
	print_cutD(cd,"temporary_displayed_cut.ps");
	system("gv temporary_displayed_cut.ps");
	system("rm temporary_displayed_cut.ps");

}

template void print_tree_graph<CutType>(CutType*, char const*);
template void print_tree_graph(Cut_Part_D_Dims*, char const*);

}

