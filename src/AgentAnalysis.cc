#include "Agent.hh"
#include "stdio.h"
#include <stdlib.h> 
#include <png.h>
#include <zlib.h>
#include <string>

void Agent::CreateAgentFromFile(char *despath,int agentidnr)
{
  InitAgent(agentidnr,0,0);
  G->CreateGenomeFromFile(despath,agentidnr);
  N->BuildNetwork(G);
  DevelopAgent();
  DetermineFitness(0);
}

void Agent::InitAgentFromFile(char *despath,int agentidnr)
{
  InitAgent(agentidnr,0,0);
  G->CreateGenomeFromFile(despath,agentidnr);
}

void Agent::FinishAgentFromFile()
{
  N->BuildNetwork(G);
  DevelopAgent();
  DetermineFitness(0);
}


void Agent::DetermineGenotypeSimilarity(Agent *M,int *genosim)
{
  Genome::iter it1;
  Genome::iter it2;
  *genosim=0;

  it1=G->ChromBBList->begin();
  it2=M->G->ChromBBList->begin();

  if(G->ChromBBList->size() != M->G->ChromBBList->size())
    *genosim=1;
  else
    {
      while(it1!=G->ChromBBList->end() || it2!=M->G->ChromBBList->end())
	{
	  if((**it1)!=(**it2))
	    {
	      *genosim=1;
	      break;
	    }
	    
	  it1++;
	  it2++;
	}
    }
}
void Agent::CloneAgentFromAgentForSafety(Agent *A,int agentidnr,int flag)
{
  InitAgent(agentidnr,0,0);
  G->CloneGenome(A->G);
  N->BuildNetwork(G);
}
/*
int Agent::DifferentExpressionPattern(double p1[NrGeneTypes],double p2[NrGeneTypes])
{
  int different;
  int i;

  //0: the same
  //1: different
  different=0;
  for(i=NrMatGeneTypes;i<NrGeneTypes;i++)
    {
      if(G->idimpact[i]==1)//this gene type/protein contributes to cell id
	{
	  if( (p1[i]<ThOff && p2[i]>=ThOn) || (p1[i]>=ThOn && p2[i]<ThOff) )
	    {
	      different=1;
	      break;
	    }
	}
    }
  return different;
}*/

void Agent::FindProteinOscillations(int *array)
{
  int counter=0, percounter=0;
  int period=0;
  int flag=0;
  int min=100, max=0;
  
  for(int i=0; i<NrGeneTypes; i++)
  {
    array[i]=0;
    max=0;
    min=100;
    for(int j=NrDevSteps/StorageInt-50; j>=50; j--)
    {
      //printf(" gene %d, E is %d\n", i, E[j][0][i]);
      if(E[j][1][i]>max){
	max=E[j][1][i];
      }
      else if (E[j][1][i]<min && E[j][1][i]>0){
	min=E[j][1][i];
      }
    }
    
    if (max-min>15 && min<100 && max>0)
      array[i]=1;
        
    
  }
  
  
}

int Agent::DifferentEquilibrium(double p1[NrGeneTypes],double p2[NrGeneTypes])
{
  int different;
  int difference;
  int d;
  int i;


  different=0;
  difference=0;
  for(i=NrMatGeneTypes;i<NrGeneTypes;i++)
    {
      d=abs(p1[i]-p2[i]);
      if(d>30)
	{
	  difference=1000;
	  break;
	}
      else
	difference+=d;
    }
  if(difference>60)
    different=1;

  return different;
}

int Agent::CompareAgents(Agent *A,int c)
{
  int i,j;
  int counter=0;
  list<Cell>::iterator iter, Aiter;

  if(c==1)//compare exact expression of all genes
    {
      for(iter=cells.begin(), Aiter=A->cells.begin();iter!=cells.end();++iter, ++Aiter)  
	  for(j=0;j<NrGeneTypes;j++)
	  {
	     if( ( (*iter).proteinstates[j]<ThOff && (*Aiter).proteinstates[j]>ThOff ) ||
	       ( (*iter).proteinstates[j]>ThOn && (*Aiter).proteinstates[j]<ThOn ) ||
	       ( ( (*iter).proteinstates[j]>ThOff && (*iter).proteinstates[j]<ThOn)  &&
	       ((*Aiter).proteinstates[j]<ThOff || (*Aiter).proteinstates[j]>ThOn)))
	      counter++;
	  }
	  
    }
    
    
  else if(c==2)//compare exact expression of segm and id genes only
    {
     printf("invalid C option: not checking domains\n");
     return 1;
    }
  else if(c==3)//compare nr of segm and domains
    {
      printf("invalid C option: not checking domains\n");
      return 1;
    }
  else if(c==4)//compare nr of segm
    {
#ifdef HOMOGENEOUSSEGMENTS
      if(nrhomogbands!=A->nrhomogbands)
#else
      if(nrlongbands!=A->nrlongbands)	
#endif
      	counter+=6;
    }
  else if(c==5)//compare nr of dom
    {
      printf("invalid C option: not checking domains\n");
      return 1;
    }
  else if(c==6)//compare whether both have segm
    {
#ifdef HOMOGENEOUSSEGMENTS
      if(A->nrhomogbands<=4 && nrhomogbands!=A->nrhomogbands)
	counter+=6;
      else if(A->nrhomogbands==5 && nrhomogbands<4)
	counter+=6;
      else if(A->nrhomogbands>5 && nrhomogbands<A->nrhomogbands-2)	      
	counter+=6;
#else
      if(A->nrlongbands<=4 && nrlongbands!=A->nrlongbands)
	counter+=6;
      else if(A->nrlongbands==5 && nrlongbands<4)
	counter+=6;
      else if(A->nrlongbands>5 && nrlongbands<A->nrlongbands-2)	      
	counter+=6;      
#endif
	      /*
#ifdef HOMOGENEOUSSEGMENTS
      if(A->nrhomogbands==2 && nrhomogbands<2)
#else
      if(A->nrlongbands==2 && nrlongbands<2)
#endif
	counter+=6;
#ifdef HOMOGENEOUSSEGMENTS
      else if(A->nrhomogbands==3 && nrhomogbands<3)
#else
      else if(A->nrlongbands==3 && nrlongbands<3)
#endif
	counter+=6;
#ifdef HOMOGENEOUSSEGMENTS
      else if(A->nrhomogbands>=4 && (nrhomogbands<(int)((double)A->nrhomogbands/2.0+0.5) || nrhomogbands<3))
#else
      else if(A->nrlongbands>=4 && (nrlongbands<(int)((double)A->nrlongbands/2.0+0.5) || nrlongbands<3))
#endif
	counter+=6;
	      */
    }
  else if(c==7)//compare whether both have same dom (once segm gene removed!!!!)
    {
      printf("invalid C option: not checking domains\n");
      return 1;
    }
  else if(c==8)//compare exact expression pattern of segm genes
    {
      //compare expression of segmentation determining gene
      j=SegmGeneNr;
      for(iter=cells.begin(), Aiter=A->cells.begin();iter!=cells.end();++iter, ++Aiter)  
	{
	     if(( (*iter).proteinstates[j]<ThOff && (*Aiter).proteinstates[j]>ThOff) ||
	       ( (*iter).proteinstates[j]>ThOn && (*Aiter).proteinstates[j]<ThOn) ||
	       ( ( (*iter).proteinstates[j]>ThOff && (*iter).proteinstates[j]<ThOn)  &&
	       ( (*Aiter).proteinstates[j]<ThOff || (*Aiter).proteinstates[j]>ThOn) ) )
	     {
	       counter++;
	       //printf("difference!\n");
	     }
	}
	
    }
    
  else if(c==9)//compare exact expression pattern of id genes
    {
     printf("invalid C option: not checking domains\n");
     return 1;
    }
 

  if(counter>5)
    return 1;
  else
    return 0;
}

void Agent::WriteGenome(int c)
{
  FILE *f;
  char fname[800];
  Genome::iter it;
  Gene *gene;
  TFBS *tfbs;
  int teller;
  double YY;
  double step=0.75;

  
  if(c==0)
    sprintf(fname,"%s/CodedGenomeAgent%.10d%s",despath,agentid,"original");
  else 
    sprintf(fname,"%s/CodedGenomeAgent%.10d_%d",despath,agentid,c);


  
  f=fopen(fname,"w");
  for(it=G->ChromBBList->begin();it!=G->ChromBBList->end();it++)
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  fprintf(f,"%i ",0);//to indicate that it is a gene
	  fprintf(f,"%i ",(gene->type));//type of gene
	}
      else if(G->IsTFBS(*it))
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  if(tfbs->weight<0)
	    fprintf(f,"%i ",-1);//to indicate that it is a repressing tfbs
	  else
	    fprintf(f,"%i ",1);//to indicate that it is a activating tfbs
	  fprintf(f,"%i ",(tfbs->type));//type of tfbs
	}
    }
      
  fclose(f);

  if(c==0)
    sprintf(fname,"%s/GenomeAgent%.10d%s.dot",despath,agentid,"original");
  else 
    sprintf(fname,"%s/GenomeAgent%.10d_%d.dot",despath,agentid,c);
   
  f=fopen(fname,"w");
  fprintf(f,"graph genome_%i {\n",agentid);
  fprintf(f,"size=\"10,4\";\n");
  teller=0;
  for(it=G->ChromBBList->begin();it!=G->ChromBBList->end();it++)
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  YY=step*teller;
	  fprintf(f,"\"%i\" [label=\"%i\",shape=rectangle,color=blue,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*gene).type,YY,0.);
	}
      else if(G->IsTFBS(*it))
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  YY=step*teller;
	  if((*tfbs).type<NrMatGeneTypes)
	    {
	      if((*tfbs).weight==1)
		fprintf(f,"\"%i\" [label=\"%i\",shape=octagon,color=green,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	      else
		fprintf(f,"\"%i\" [label=\"%i\",shape=octagon,color=red,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	    }
	  else if((*tfbs).type<(NrMatGeneTypes+NrSignGeneTypes))
	    {
	      if((*tfbs).weight==1)
		fprintf(f,"\"%i\" [label=\"%i\",shape=diamond,color=green,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	      else
		fprintf(f,"\"%i\" [label=\"%i\",shape=diamond,color=red,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	    }
	  else
	    {
	      if((*tfbs).weight==1)
		fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=green,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	      else
		fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=red,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	    }
	    
	}
      teller++;
    }
  fprintf(f,"}\n");
  fclose(f);
}

int mapEE(double EE)
{
  int val;
  
  
  
  return val;
}
static string NumToHex(int val)
{
  char hex[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
  char hold[]={'0','0','0'};
  int rep=0;
  int num=val;
  int rest=num%16;
  char reverse[]={'0','0','0'};
  
  hold[rep]=hex[rest]; 
  while((num-rest)!=0){
    rep++;
    num=(num-rest)/16;
    rest=num%16;
    hold[rep]=hex[rest]; 
  }
  
  reverse[0]=hold[1];
  reverse[1]=hold[0];
  reverse[2]='\0';
  
  return string(reverse);
}

static string ExprToHexColor(double val, double minval, double maxval, int code)
{
  int R,G,B;
  string RGB;
  double minsize, maxsize,nval;
  ostringstream strs;
  
  double conval=((val-minval)/(maxval-minval))*255.;
  switch (code){
    case 1: //expression
      if(conval<maxval){
	R=(int)max(255.-conval*0.5,0.);//R
	G=(int)max(255.-conval*1.5,0.);//G
	B=(int)max(255.-conval*3.,0.);//B 
      }
      else{
	R=125;
	G=0;
	B=0;
      }
      RGB=NumToHex(R)+NumToHex(G)+NumToHex(B);
      break;
    case 2: //DD value of a gene to color
      R=(int)max(255.-conval,0.);
      G=(int)max(255.-conval,0.);
      B=255;
      RGB=NumToHex(R)+NumToHex(G)+NumToHex(B);
      break;
    case 3: //EE value of a gene to font size
      minsize=8.; 
      maxsize=24.;
      conval/=255.;
      nval=minsize+conval*(maxsize-minsize);
      strs <<nval;
      RGB=strs.str();
      break;
    case 4://HH value of a TFBS to an edge thickness
      minsize=0.5;
      maxsize=5.;
      conval/=255.;
      nval=minsize+conval*(maxsize-minsize);
      strs <<nval;
      RGB=strs.str();
      break;
    default:
      printf("ExprToHexColor: error: wrong code\n");
      break;
      
  }
  return RGB;
}

void Agent::WriteNetwork(char *c)
{
  FILE *f;
  char fname[800];
  list<Vertex *> *VL;
  list<list <Edge *> > *AL;
  Edge *E;
  Vertex *V;
  Network::iterv iv;
  Network::iterel el;
  Network::itere e;
  int nodelabelin;
  int nodelabelout;
  Gene *gene;
  int weight;
  char ss[200];
  string ee;
  string dd;
  string hh;
  VL=(*N).VL;
  AL=(*N).AL;
 
  sprintf(fname,"%s/NetworkAgent%.10d_%s.dot",despath,agentid,c);
  
  f=fopen(fname,"w");
  fprintf(f,"digraph network_%i {\n",agentid);
  fprintf(f,"size=\"10,10\";\n");
  fprintf(f,"overlap=false; \n");
  fprintf(f,"splines=true; \n");
  
  iv=VL->begin();
  for(el=AL->begin();el!=AL->end();el++)
    {
      V=(*iv);
      nodelabelin=(*V).label;
      gene=(*V).Gen;
      if(gene->type==0)
	{
	  sprintf(ss,"5,0");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==1)
	{ 
	  sprintf(ss,"1,2");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==2)
	{ 
	  sprintf(ss,"3,2");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==3)
	{ 
	  sprintf(ss,"5,2");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==4)
	{ 
	  sprintf(ss,"7,2");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
}
      else if(gene->type==5)
	{ 
	  sprintf(ss,"9,2");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==6)
	{ 
	  sprintf(ss,"1,4");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==7)
	{ 
	  sprintf(ss,"3,4");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==8)
	{ 
	  sprintf(ss,"5,4");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==9)
	{ 
	  sprintf(ss,"7,4");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==10)
	{ 
	  sprintf(ss,"9,4");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==11)
	{ 
	  sprintf(ss,"1,6");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==12)
	{ 
	  sprintf(ss,"3,6");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==13)
	{ 
	  sprintf(ss,"5,6");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==14)
	{ 
	  sprintf(ss,"7,6");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else if(gene->type==15)
	{ 
	  sprintf(ss,"9,6");
	  ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	  dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	  fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\", pos=\"%s\",pin=true];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str(), ss);
	}
      else
      {
	ee=ExprToHexColor(gene->EE, Emin, Emax,3);
	dd=ExprToHexColor(gene->DD, Dmin, Dmax,2);
	fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,fontsize=%s, style=filled, fillcolor=\"#%s\"];\n",nodelabelin,(*gene).type,ee.c_str(),dd.c_str());
      }
      for(e=el->begin();e!=el->end();e++)
	{
	  E=(*e);
	  V=(*E).V;
	  nodelabelout=(*V).label;
	  weight=(*E).weight;
	  hh=ExprToHexColor((*E).HH,Hmin,Hmax, 4);
	  if(weight==1)
	    {
	      fprintf(f,"\"%i\" -> \"%i\" [color=green,penwidth=%s];\n",nodelabelout,nodelabelin,hh.c_str());
	      
	    }
	  else if(weight==-1)
	    {
	      fprintf(f,"\"%i\" -> \"%i\" [color=red,penwidth=%s];\n",nodelabelout,nodelabelin,hh.c_str());	
	    }
	}
      iv++;
    }
  fprintf(f,"}\n");
  fclose(f);
}




void Agent::WriteNetworkProgression(int c, char* destdir, int cell)
{
  FILE *f;
  char fname[800], fname2[800];
  list<Vertex *> *VL;
  list<list <Edge *> > *AL;
  Edge *E;
  Vertex *V;
  Network::iterv iv;
  Network::iterel el;
  Network::itere e;
  int nodelabelin;
  int nodelabelout;
  Gene *gene;
  int weight;
  char ss[200];
  string expr;

  VL=(*N).VL;
  AL=(*N).AL;
 
  if(c==0)
    sprintf(fname,"%s/NetworkCell%d%s_",destdir,cell,"original");
  else if(c==1)
    sprintf(fname,"%s/NetworkCell%d%s_",destdir,cell, "core");
  else if(c==2)
   {
    printf("Wrong option. no domains\n");
    return;
   }
  else if(c==3)
   {
    printf("Wrong option. no domains\n");
    return;
   }
  else if(c==4)
    sprintf(fname,"%s/NetworkCell%d%s_",destdir,cell,"nrsegm");
  else if(c==5)
    {
    printf("Wrong option. no domains\n");
    return;
    }
  else if(c==6)
    sprintf(fname,"%s/NetworkCell%d%s_",destdir,cell,"minsegm");
  else if(c==7)
   {
    printf("Wrong option. no domains\n");
    return;
   }
  else if(c==8)
    sprintf(fname,"%s/NetworkCell%d%s_",destdir,cell,"segm");
  else if(c==9)
   {
    printf("Wrong option. no domains\n");
    return;
   }
  else if(c==10)
    sprintf(fname,"%s/NetworkCell%d%s_",destdir,cell,"coreminsegm");
  
  char number[50];
  int xx,yy;
  for (int i=0; i<NrStorages;i++){
    sprintf(number,"%04d.dot",i);
    strcpy(fname2,fname);
    strcat(fname2, number);
    f=fopen(fname2,"w");
    fprintf(f,"digraph network_%i {\n",agentid);
    fprintf(f,"size=\"10,10\";\n");
    fprintf(f,"overlap=false; \n");
    fprintf(f,"splines=true; \n");
    
    
    iv=VL->begin();
    for(el=AL->begin();el!=AL->end();el++)
    {
      V=(*iv);
      nodelabelin=(*V).label;
      gene=(*V).Gen;
      if(gene->type==0)
      {
	expr=ExprToHexColor((double)cellpattern[i][cell/10][gene->type],0.,Emax/Dmin,1);
	sprintf(ss,"5,0");
	fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,pos=\"%s\",pin=true, style=filled, fillcolor=\"#%s\"];\n",nodelabelin,(*gene).type,ss,expr.c_str());
      }
      else if(gene->type<=15)
      {
	xx=((gene->type-1)%5)+1;
	yy=(gene->type/5)*2+2;
	sprintf(ss,"%d,%d",xx+xx-1,yy);
	expr=ExprToHexColor((double)cellpattern[i][cell/10][gene->type],0.,Emax/Dmin,1);
	fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold,pos=\"%s\",pin=true, style=filled, fillcolor=\"#%s\"];\n",nodelabelin,(*gene).type,ss,expr.c_str());
      }      
      else
	fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=black,style=bold];\n",nodelabelin,(*gene).type);
      
      for(e=el->begin();e!=el->end();e++)
      {
	E=(*e);
	V=(*E).V;
	nodelabelout=(*V).label;
	weight=(*E).weight;
	if(weight==1)
	{
	  fprintf(f,"\"%i\" -> \"%i\" [color=green,penwidth=3];\n",nodelabelout,nodelabelin);
	  
	}
	else if(weight==-1)
	{
	  fprintf(f,"\"%i\" -> \"%i\" [color=red,penwidth=3];\n",nodelabelout,nodelabelin);	
	}
      }
      iv++;
    }
    fprintf(f,"}\n");
    fclose(f);
  }
}


/*
void Agent::WriteMovieEmbryology()
{
  int i,j,k,l;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[800];
  const int L1=NrFinalCells;
  const int L2=10;
  const int Lt=4*L1+3*L2;
  const int W1=20;
  const int W2=6;//5
  const int Wt=4*W1+3*W2;
  int schema[Lt][Wt];
  unsigned char RGBdata[Lt*Wt*3];
  char s[800];

  if(filecounter==0)
    {
      //make dir
      sprintf(s,"mkdir %s/MovieAgent%.10d",despath,agentid);
      system(s);
      printf("made directory\n");
    }
  sprintf(fname,"%s/MovieAgent%.10d/snapshot%.10d.png",despath,agentid,filecounter+20);

  
  for(i=0;i<Lt;i++)
    for(j=0;j<Wt;j++)for(i=0;i<NrFinalCells;i++)
	{
	  if((CL[i]->proteinstates[j]<ThOff && A->CL[i]->proteinstates[j]>ThOff) ||
	     (CL[i]->proteinstates[j]>ThOn && A->CL[i]->proteinstates[j]<ThOn) ||
	     ((CL[i]->proteinstates[j]>ThOff && CL[i]->proteinstates[j]<ThOn)) &&
	     ((A->CL[i]->proteinstates[j]<ThOff || A->CL[i]->proteinstates[j]>ThOn)))
	    counter++;for(i=0;i<NrFinalCells;i++)
	{
	  if((CL[i]->proteinstates[j]<ThOff && A->CL[i]->proteinstates[j]>ThOff) ||
	     (CL[i]->proteinstates[j]>ThOn && A->CL[i]->proteinstates[j]<ThOn) ||
	     ((CL[i]->proteinstates[j]>ThOff && CL[i]->proteinstates[j]<ThOn)) &&
	     ((A->CL[i]->proteinstates[j]<ThOff || A->CL[i]->proteinstates[j]>ThOn)))
	    counter++;
	}
	}
      schema[i][j]=0; for(i=0;i<NrFinalCells;i++)
	{
	  if((CL[i]->proteinstates[j]<ThOff && A->CL[i]->proteinstates[j]>ThOff) ||
	     (CL[i]->proteinstates[j]>ThOn && A->CL[i]->proteinstates[j]<ThOn) ||
	     ((CL[i]->proteinstates[j]>ThOff && CL[i]->proteinstates[j]<ThOn)) &&
	     ((A->CL[i]->proteinstates[j]<ThOff || A->CL[i]->proteinstates[j]>ThOn)))
	    counter++;
	}

  //Settings for Seed7:segments first:
  
  /*
  //gene 0, maternal gene expression
  for(i=0;i<L1;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=CL[i]->proteinstates[12];
  //gene 1
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=CL[i-(L1+L2)]->proteinstates[13];
  //gene 2
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=CL[i-(2*L1+2*L2)]->proteinstates[14];
  //gene 3
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=2*CL[i-(3*L1+3*L2)]->proteinstates[15];

  //gene 4
  for(i=0;i<L1;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=1.5*CL[i]->proteinstates[8];
  //gene 5, segmentation gene
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=1.5*CL[i-(L1+L2)]->proteinstates[9];
  //gene 6
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=1.5*CL[i-(2*L1+2*L2)]->proteinstates[10];
  //gene 7
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=1.5*CL[i-(3*L1+3*L2)]->proteinstates[11];
 
  
  //gene 8
  for(i=0;i<L1;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=CL[i]->proteinstates[4];
  //gene 9
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=2*CL[i-(L1+L2)]->proteinstates[5];
  //gene 10
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=3.5*CL[i-(2*L1+2*L2)]->proteinstates[6];
  //gene 11
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=2*CL[i-(3*L1+3*L2)]->proteinstates[7];

  
  //gene 12
  for(i=0;i<L1;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=2.5*CL[i]->proteinstates[0];
  //gene 13
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=2.5*CL[i-(L1+L2)]->proteinstates[1];
  //gene 14
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=4.5*CL[i-(2*L1+2*L2)]->proteinstates[2];
  //gene 15
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=2.5*CL[i-(3*L1+3*L2)]->proteinstates[3];
 
  */
  /*
  list<Cell>::iterator iter;
  //New settings for seed 3, segments simultaneous
 
  //gene 0, maternal gene expression
  for(i=0;i<L1;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i]->proteinstates[12]));
  //gene 1
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(L1+L2)]->proteinstates[13]));
  //gene 2
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(2*L1+2*L2)]->proteinstates[14]));
  //gene 3
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=0;j<W1;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(3*L1+3*L2)]->proteinstates[15]));

  //gene 4
  for(i=0;i<L1;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i]->proteinstates[8]));
  //gene 5, segmentation gene
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(L1+L2)]->proteinstates[9]));
  //gene 6
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(2*L1+2*L2)]->proteinstates[10]));
  //gene 7
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=W1+W2;j<2*W1+W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(3*L1+3*L2)]->proteinstates[11]));
 
  
  //gene 8
  for(i=0;i<L1;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i]->proteinstates[4]));
  //gene 9
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=minint(255,(int)(4*CL[i-(L1+L2)]->proteinstates[5]));
  //gene 10
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(2*L1+2*L2)]->proteinstates[6]));
  //gene 11
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=2*W1+2*W2;j<3*W1+2*W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(3*L1+3*L2)]->proteinstates[7]));

  
  //gene 12
  for(i=0;i<L1;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i]->proteinstates[0]));
  //gene 13
  for(i=L1+L2;i<2*L1+L2;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(L1+L2)]->proteinstates[1]));
  //gene 14
  for(i=2*L1+2*L2;i<3*L1+2*L2;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(2*L1+2*L2)]->proteinstates[2]));
  //gene 15
  for(i=3*L1+3*L2;i<4*L1+3*L2;i++)
    for(j=3*W1+3*W2;j<4*W1+3*W2;j++)
      schema[i][j]=minint(255,(int)(2.5*CL[i-(3*L1+3*L2)]->proteinstates[3]));
 





  for(i=0;i<Lt;i++)
    for(j=0;j<Wt;j++)
      {
	RGBdata[j*Lt*3+i*3+0]=schema[i][j];
	RGBdata[j*Lt*3+i*3+1]=schema[i][j];
	RGBdata[j*Lt*3+i*3+2]=schema[i][j];
      }


  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,Lt,Wt,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = Wt-1; row >= 0; row-- ) 
    {
      row_pointer = ( RGBdata + Lt * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);

  filecounter++;

  printf("wrote png file %i of agent %i\n",filecounter,agentid);
  
}
*/
static int AgeToColor(int minage,int maxage, int age)
{
  int color;
  
  if(age<0)
    return -1;
  
  color=(int)((double)(age-minage)/(double)(maxage-minage)*255.);
 
  if(color>255)
    color=255;
  
  
  return color;
}
  
void Agent::WriteFullAgeProfile(char *dirname,int c)
{
  int i,j;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*NrStorages;
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  // int type;
  int color;
  
  
    for(i=0;i<NrStorages;i++)
  {
    for(j=0;j<NrFinalCells;j++)
    {
      color=AgeToColor(0,NrDevSteps,ages[i][j]);
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<zoom;jj++)
	  celltypes[zoom*i+ii][zoom*j+jj]=color;
    }
  }

  for(i=0;i<LL;i++)
    for(j=0;j<WW;j++)
    {
      if(celltypes[i][j]>=0) //graded coloring of cell age
      {
	RGBdata[i*WW*3+j*3+0]=max(0.0, 255.-0.5*celltypes[i][j]);
	RGBdata[i*WW*3+j*3+1]=max(0.0, 255.-celltypes[i][j]);
	RGBdata[i*WW*3+j*3+2]=max(0.0, 255.-3*celltypes[i][j]);
	//gnuplot2 gradient; makes no difference. is worse if anything.
// 	RGBdata[i*WW*3+j*3+0]=max(min(255.,(double)(celltypes[i][j]-64)*4.),0.);//R	
// 	RGBdata[i*WW*3+j*3+1]=max(min(255.,(double)(celltypes[i][j]-128)*4.),0.);//G
// 	if(celltypes[i][j]<128) //blue is more complicated
// 	  RGBdata[i*WW*3+j*3+2]=max(min(255.,(double)celltypes[i][j]*4.),0.);
// 	else if(celltypes[i][j]<192)
// 	  RGBdata[i*WW*3+j*3+2]=max(min(255.,255.-(double)(celltypes[i][j]-127)*4.),0.);
// 	else
// 	  RGBdata[i*WW*3+j*3+2]=max(min(255.,(double)(celltypes[i][j]-192)*4.),0.);
      }
      
      
      else //non-tissue part
      {
	RGBdata[i*WW*3+j*3+0]=127;
	RGBdata[i*WW*3+j*3+1]=127;
	RGBdata[i*WW*3+j*3+2]=127;
	
      }
    }

  sprintf(fname,"%s/%s/Age%.10d_%d.png",despath,dirname,agentid, c);
  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,WW,LL,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
//   // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);
}

void Agent::WriteTimepointAgeProfile(char *dirname, int c, int timepoint)
{
  int i,j;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*4;
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  // int type;
  int color;
  int minage=1000;
  int maxage=0;
  
 // printf("timepoint %d, ages:\n", timepoint);
  
  for(j=0;j<NrFinalCells;j++)
  {
    //printf("%d ",ages[timepoint][j]);
    if(ages[timepoint][j]<minage && ages[timepoint][j]>=0)
      minage=ages[timepoint][j];
    if(ages[timepoint][j]>maxage)
      maxage=ages[timepoint][j];
  }
 // printf("\n\n");
 // printf("minage: %d, maxage=%d\n",minage, maxage);
   
  
  for(i=0;i<4;i++)
  {
    for(j=0;j<NrFinalCells;j++)
    {
      color=AgeToColor(minage,maxage, ages[timepoint][j]);
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<zoom;jj++)
	  celltypes[zoom*i+ii][zoom*j+jj]=color;
    }
  }

  for(i=0;i<LL;i++)
    for(j=0;j<WW;j++)
    {
      if(celltypes[i][j]>=0) //graded coloring of cell age
      {
	RGBdata[i*WW*3+j*3+0]=max(0.0, 255.-0.5*celltypes[i][j]);
	RGBdata[i*WW*3+j*3+1]=max(0.0, 255.-celltypes[i][j]);
	RGBdata[i*WW*3+j*3+2]=max(0.0, 255.-3*celltypes[i][j]);
// 	RGBdata[i*WW*3+j*3+0]=max(min(255.,(double)(celltypes[i][j]-64)*4.),0.);//R	
// 	RGBdata[i*WW*3+j*3+1]=max(min(255.,(double)(celltypes[i][j]-128)*4.),0.);//G
// 	if(celltypes[i][j]<128) //blue is more complicated
// 	  RGBdata[i*WW*3+j*3+2]=max(min(255.,(double)celltypes[i][j]*4.),0.);
// 	else if(celltypes[i][j]<192)
// 	  RGBdata[i*WW*3+j*3+2]=max(min(255.,255.-(double)(celltypes[i][j]-127)*4.),0.);
// 	else
// 	  RGBdata[i*WW*3+j*3+2]=max(min(255.,(double)(celltypes[i][j]-192)*4.),0.);
      }
      
      else //non-tissue part
      {
	RGBdata[i*WW*3+j*3+0]=127;
	RGBdata[i*WW*3+j*3+1]=127;
	RGBdata[i*WW*3+j*3+2]=127;
	
      }
    }

  sprintf(fname,"%s/%s/AgeAt%d_%.10d_%d.png",despath,dirname,timepoint, agentid, c);
  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,WW,LL,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
//   // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);
}

static int typecolormap[18][3]=
{  {255, 0, 0},//red 0
   {0, 255, 0},//green 1 
   {0, 0, 255},//blue 2
   {255, 153, 51}, //orange 3(not the 'perfect' orange but 1 shade lighter for better distinction)
   {0, 255, 255},//cyan 4
   {255, 0, 255},//magenta 5
   {255, 255, 0}, //yellow 6
   {153, 51, 255}, //purple 7 (same as orange)
   {205,92,92}, //indian red 8
   {0, 128, 0}, //darkgreen 9
   {218,165,32}, //goldenrod 10
   {65,105,225}, //royalblue, 11
   {139,0,0}, //dark red 12
   {102,205,170}, //medium aquamarine 13
   {210,105,30}, //chocolate 14
   {147,112,219}, //light purple 15
   {255, 255, 255}, //white 16
   {0, 0, 0}//black 17
  };
  
void Agent::WriteCloneProfile(char *dirname, int c, int timepoint)
{
  int i,j;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*NrStorages;
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  // int type;
  int color;
  int nrcells1=0, nrcells2=0;
  
  //the crucial time points are the second and third; min age at any consecutive time point is at least timepoint-2 or

 int counter;
  for(i=0;i<NrStorages;i++)
  {
    counter=-1;
    for(j=0;j<NrFinalCells;j++)
    {
      if(ages[i][j]<0) //not a cell
	color=-1;
      else if(i<timepoint)//before the interesting point, all cells are black
	color=17;
      else if((ages[i][j]/StorageInt)>=i-timepoint && j>=InitNrCells-1)//switch colour if we have an original cell, don't count head.
      {
	counter++;
	color=counter%18;
      }
      else
	color=counter%18;
	
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<zoom;jj++)
	  celltypes[zoom*i+ii][zoom*j+jj]=color;
    }
  }

  for(i=0;i<LL;i++)
    for(j=0;j<WW;j++)
    {
      if(celltypes[i][j]>=0) //graded coloring of cell age
      {
	RGBdata[i*WW*3+j*3+0]=typecolormap[celltypes[i][j]][0];
	RGBdata[i*WW*3+j*3+1]=typecolormap[celltypes[i][j]][1];
	RGBdata[i*WW*3+j*3+2]=typecolormap[celltypes[i][j]][2];
      }
      
      else //non-tissue part
      {
	RGBdata[i*WW*3+j*3+0]=127;
	RGBdata[i*WW*3+j*3+1]=127;
	RGBdata[i*WW*3+j*3+2]=127;
	
      }
    }

  sprintf(fname,"%s/%s/ClonesPoint%d_%.10d_%d.png",despath,dirname,timepoint,agentid, c);
  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,WW,LL,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
//   // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);
}

void Agent::WriteCloneTimecourse(char *subdir, int c, int timepoint)
{
  //file opening 
  FILE *f;
  char fname[500];
  sprintf(fname,"%s/%s/Clonetimecourse%d_%d",despath, subdir,timepoint, c);
  f=fopen(fname, "w");
  
  for(int k=0; k<NrGeneTypes; k++)
  {
    for(int i=timepoint; i<NrStorages; i++)
    {
      fprintf(f,"%d\t",i);
      for(int j=0; j<NrFinalCells; j++)
      {
	if((ages[i][j]/StorageInt)>=i-timepoint && j>=InitNrCells-1 && ages[i][j]>=0)
	  fprintf(f, "%d\t",E[i][j][k]);
      }
      fprintf(f,"\n");
    }
    fprintf(f,"\n\n\n");
  }
  fclose(f);
}

void Agent::WriteGeneEmbryology(char *subdir, int c, int gene)
{
  int i,j,k,l;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*(NrStorages);
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  int type;
  int color;

  const int seg=gene;
  int max=0;
  //find max expression of this gene
  for(i=0;i<NrStorages;i++)
    for(j=0;j<NrFinalCells;j++)
      if(E[i][j][seg]>max)
	max=E[i][j][seg];
  
  for(i=0;i<NrStorages;i++)
    for(j=0;j<NrFinalCells;j++)
    {
      color=SegmentToColor(E[i][j][seg],max);
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<zoom;jj++)
	  celltypes[zoom*i+ii][zoom*j+jj]=color;
    }

    //blue to white scale
    for(i=0;i<LL;i++)
      for(j=0;j<WW;j++)
      {
	RGBdata[i*WW*3+j*3+0]=255-celltypes[i][j];
	RGBdata[i*WW*3+j*3+1]=255-celltypes[i][j];
	RGBdata[i*WW*3+j*3+2]=255;
      }

 
  if(c==0)
    sprintf(fname,"%s/%s/Gene%dAgent%.10d%s.png",despath,subdir,gene, agentid,"original");
  else
    sprintf(fname,"%s/%s/Gene%dAgent%.10d_%d.png",despath,subdir,gene, agentid,c);
    

  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,WW,LL,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);

  //printf("spacetimeplot development png written to file\n");

}

void Agent::WriteTimepointGene(char *subdir, int c, int gene, int timepoint)
{
  int i,j,k,l;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*4;
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  int type;
  int color;

  int max=0;
  //find max expression of this gene
  for(j=0;j<NrFinalCells;j++)
    if(E[timepoint][j][gene]>max)
      max=E[timepoint][j][gene];
  
  for(i=0;i<4;i++)
    for(j=0;j<NrFinalCells;j++)
    {
      color=SegmentToColor(E[timepoint][j][gene],max);
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<zoom;jj++)
	  celltypes[zoom*i+ii][zoom*j+jj]=color;
    }

    //blue to white scale
    for(i=0;i<LL;i++)
      for(j=0;j<WW;j++)
      {
	RGBdata[i*WW*3+j*3+0]=255-celltypes[i][j];
	RGBdata[i*WW*3+j*3+1]=255-celltypes[i][j];
	RGBdata[i*WW*3+j*3+2]=255;
      }

 
  if(c==0)
    sprintf(fname,"%s/%s/Gene%dAt%d_Agent%.10d%s.png",despath,subdir,gene,timepoint, agentid,"original");
  else
    sprintf(fname,"%s/%s/Gene%dAt%d_Agent%.10d_%d.png",despath,subdir,gene,timepoint, agentid,c);
    

  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,WW,LL,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);

  //printf("spacetimeplot development png written to file\n");

}
void Agent::WriteVideo(int gene)
{
  int genemax=0, segmax=0;
  //find max expression for scaling
  for (int i=0; i<NrStorages; i++)
    for (int j=0; j<NrFinalCells; j++)
    {
      if(genemax<E[i][j][gene])
	genemax=E[i][j][gene];
      if(segmax<E[i][j][SegmGeneNr])
	segmax=E[i][j][SegmGeneNr];
    }
  for(int i=0; i<NrStorages; i++)
  {
    WriteVideoPicture(gene, i, genemax);
    WriteVideoPicture(SegmGeneNr, i, segmax);
  }
}

void Agent::WriteVideoPicture(int gene, int time, int max)
{
  int i,j,k,l;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int width=zoom*3;
  const int WW=zoom*NrFinalCells;
  int celltypes[width][WW];
  unsigned char RGBdata[width*WW*3];
  int type;
  int color;

  const int seg=gene; //the gene to make a picture of.
  //printf("am here!!\n");
  int tissuelength=0;

  if(cells.size()<NrFinalCells)
    tissuelength=cells.size()*zoom;
  else
    tissuelength=-1;
  
  for(j=0;j<NrFinalCells;j++)
  {
    color=SegmentToColor(E[time][j][seg],max);
    for(ii=0;ii<width;ii++)
    {
      for(jj=0;jj<zoom;jj++)
      {
	celltypes[ii][zoom*j+jj]=color;
      }
    }
  }
  if(tissuelength>0)
  {
    for(ii=0;ii<width;ii++)
    {
      celltypes[ii][tissuelength]=-1;
    }
  }
  
  //blue to white scale
  for(i=0;i<width;i++)
  {
    for(j=0;j<WW;j++)
    {
      if(celltypes[i][j]>=0)
      {
	RGBdata[i*WW*3+j*3+0]=255-celltypes[i][j];
	RGBdata[i*WW*3+j*3+1]=255-celltypes[i][j];
	RGBdata[i*WW*3+j*3+2]=255;
      }
      else
	{
	RGBdata[i*WW*3+j*3+0]=0;
	RGBdata[i*WW*3+j*3+1]=0;
	RGBdata[i*WW*3+j*3+2]=0;
      }
    }
  }
  
  sprintf(fname,"%s/Gene%dVid%.10d.png",despath,gene,time);
  
  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,WW,width,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <width; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);

  //printf("spacetimeplot development png written to file\n");

}


void Agent::WriteExtraSignalProfiles(int c)
{
  
  FILE *f;
  char fname[800];
  int i,j,k;

  char ss[200];
  if(c==0)
    sprintf(ss,"original");
  else
    sprintf(ss,"_%d",c);
    
 
  for(k=0;k<NrStorages;k+=5 )//for(k=0;k<NrStorages;k+=10)
    {
      sprintf(fname,"%s/signal%.10d",despath,k);
      f=fopen(fname,"w");
      for(i=0;i<NrFinalCells;i++)
	{      
	  fprintf(f,"%i\t",i);
	  for(j=0;j<NrGeneTypes;j++)
	    fprintf(f,"%i\t",(int)E[k][i][j]);
	  fprintf(f,"\n"); 
	}
      fclose(f);
    }
    

  for(k=0;k<NrStorages;k+=10)
    {
      sprintf(fname,"%s/SignalProfiles%.10d%s_time%i",despath,agentid,ss,k*StorageInt);
      f=fopen(fname,"w");
      for(i=0;i<NrFinalCells;i++)
	{      
	  fprintf(f,"%i\t",i);
	  for(j=0;j<NrGeneTypes;j++)
	    fprintf(f,"%i\t",(int)E[k][i][j]);
	  fprintf(f,"\n"); 
	}
      fclose(f);
    }

  
 

  for(j=0;j<100;j+=5)
    {
      sprintf(fname,"%s/TemporalDynamics%.10d%s_Cell%i",despath,agentid,ss,j);
      f=fopen(fname,"w");
      for(i=0;i<NrStorages;i++)
	{
	  fprintf(f,"%i\t",i);
	  for(k=0;k<NrGeneTypes;k++)
	    fprintf(f,"%i\t",(int)E[i][j][k]);
	  fprintf(f,"\n");
	}
      fclose(f);
    }
  j=99;
  sprintf(fname,"%s/TemporalDynamics%.10d%s_Cell%i",despath,agentid,ss,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",(int)E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  
  
  for(j=0;j<NrFinalCells/10;j++)
  {
    sprintf(fname,"%s/TemporalDynamics%.10d%s_CellId%i",despath,agentid,ss,j*10);
    f=fopen(fname,"w");
    for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",cellpattern[i][j][k]);
      fprintf(f,"\n");
    }
    fclose(f);
  }
  
}


void Agent::WriteBasicProperties(int aid)
{
  FILE *f;
  char fname[800];

  sprintf(fname,"%s/%s",despath,"BasicPropertiesOfAgents");
  f=fopen(fname,"a");
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",nrbands);//2
  fprintf(f,"%i\t",nrlongbands);//3
  fprintf(f,"%i\t",nrhomogbands);//4
  fprintf(f,"%i\t",nrlonganddifbands);//5
  fprintf(f,"%i\t",nrdifbutrepeatedbands);//6
  fprintf(f,"%f\t",nonexpfitness);//7
  fprintf(f,"%f\t",fitness);//8
  fprintf(f,"%f\t",glpenalty);//9
  fprintf(f,"%f\t",instpenalty);//10
  fprintf(f,"%f\t",shortsegpenalty);//11
  fprintf(f,"%i\t",nrlongbands+nrlonganddifbands);//12
  fprintf(f,"%i\t",nrhomogbands+nrlonganddifbands);//13
  fprintf(f,"\n");
  fclose(f);
}


void Agent::DetermineGenomeAndNetworkProperties(char *filename)
{
  FILE *f;
  char fname[800];
  int i;

  Genome:: iter it;
  TFBS *tfbs;
  nrplustfbs=0;
  nrmintfbs=0;
  int indegree[200];
  int outdegree[200];
  int totdegree[200];
  int nrtfbsmatchinggenetype[NrGeneTypes];
  double indegreedistr[15];
  double outdegreedistr[15];
  double totdegreedistr[15];
  /*static double avgindegreedistr[15];
  static double avgoutdegreedistr[15];
  static double avgtotdegreedistr[15];*/ //seem to serve no purpose
  int deggenetypes[16][3]; //my addition

  for(i=0;i<200;i++)
    {
      indegree[i]=0;
      outdegree[i]=0;
      totdegree[i]=0;
     
    }
  for(i=0;i<NrGeneTypes;i++)
  {
    nrtfbsmatchinggenetype[i]=0; //outdegree
    deggenetypes[i][0]=0;//indegree (total not corrected for duplicates)
    deggenetypes[i][1]=0;//outdegree (total not corrected for duplicates)
    deggenetypes[i][2]=0;//nr of copies of this gene type
    
  }
  for(i=0;i<15;i++)
    {
      indegreedistr[i]=0;
      outdegreedistr[i]=0;
      totdegreedistr[i]=0;
    }

  int genecounter=0;  
  it=G->ChromBBList->begin();
 
  while(it!=G->ChromBBList->end())
    {
      if(G->IsTFBS(*it))//tfbs
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  if(tfbs->weight==1)
	    nrplustfbs++;
	  else
	    nrmintfbs++;

	  indegree[genecounter]++; //don't take into account that there may be multiple genes of the same type impinging on this TFBS
	  totdegree[genecounter]++;
	  nrtfbsmatchinggenetype[tfbs->type]++;
	}
      else//gene
	genecounter++;
      
      it++;
    }
 
 Gene *gene;
  genecounter=0;
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsTFBS(*it))//tfbs
	;
      else
	{
	  gene=dynamic_cast<Gene *>(*it);
	  deggenetypes[gene->type][0]+=indegree[genecounter];
	  outdegree[genecounter]=nrtfbsmatchinggenetype[gene->type];
	  deggenetypes[gene->type][1]+=outdegree[genecounter];
	  totdegree[genecounter]+=nrtfbsmatchinggenetype[gene->type];
	  deggenetypes[gene->type][2]++;
	  genecounter++;
	}
      it++;
    }
    
  //print in and out degree per gene (when duplicated, it's an average) Also note that if a gene gets input from duplicated gene, this is only counted as 1! (one TFBS)
  sprintf(fname,"%s/%s%s",despath,"DegreesPerGeneType_",filename);
  f=fopen(fname,"w");
  for(i=0; i<NrGeneTypes; i++)
    fprintf(f, "%d\t%.2lf\t%.2lf\t%d\n",i,(double)deggenetypes[i][0]/(double)deggenetypes[i][2], (double) deggenetypes[i][1]/(double)deggenetypes[i][2],deggenetypes[i][2]);
 
  avgindegree=0;
  maxindegree=0;
  minindegree=1000;
  stdindegree=0;
  avgoutdegree=0;
  maxoutdegree=0;
  minoutdegree=1000;
  stdoutdegree=0;
  avgtotdegree=0;
  maxtotdegree=0;
  mintotdegree=1000;
  stdtotdegree=0;
  for(i=0;i<genecounter;i++)
    {
      if(indegree[i]<15)
	indegreedistr[indegree[i]]++;
      else
	indegreedistr[14]++;
      if(outdegree[i]<15)
	outdegreedistr[outdegree[i]]++;
      else
	outdegreedistr[14]++;
      if(totdegree[i]<15)
	totdegreedistr[totdegree[i]]++;
      else
	totdegreedistr[14]++;

      avgindegree+=indegree[i];
      avgoutdegree+=outdegree[i];
      avgtotdegree+=totdegree[i];
      
      stdindegree+=indegree[i]*indegree[i];
      stdoutdegree+=outdegree[i]*outdegree[i];
      stdtotdegree+=totdegree[i]*totdegree[i];

      if(indegree[i]>maxindegree)
	maxindegree=indegree[i];
      if(indegree[i]<minindegree)
	minindegree=indegree[i];
      
      if(outdegree[i]>maxoutdegree)
	maxoutdegree=outdegree[i];
      if(outdegree[i]<minoutdegree)
	minoutdegree=outdegree[i];

      if(totdegree[i]>maxtotdegree)
	maxtotdegree=totdegree[i];
      if(totdegree[i]<mintotdegree)
	mintotdegree=totdegree[i];
    }
  for(i=0;i<15;i++)
    {
      indegreedistr[i]/=(double)genecounter;
      outdegreedistr[i]/=(double)genecounter;
      totdegreedistr[i]/=(double)genecounter;
    }
 /* for(i=0;i<15;i++)
    {
      avgindegreedistr[i]+=indegreedistr[i];
      avgoutdegreedistr[i]+=outdegreedistr[i];
      avgtotdegreedistr[i]+=totdegreedistr[i];
    }
    for(i=0;i<15;i++)
    {
      avgindegreedistr[i]/=100;
      avgoutdegreedistr[i]/=100;
      avgtotdegreedistr[i]/=100;
    }*/



  avgindegree/=(double)genecounter;
  avgoutdegree/=(double)genecounter;
  avgtotdegree/=(double)genecounter;

  stdindegree/=(double)genecounter;
  stdindegree=stdindegree-avgindegree*avgindegree;
  stdindegree=sqrt(stdindegree);

  stdoutdegree/=(double)genecounter;
  stdoutdegree=stdoutdegree-avgoutdegree*avgoutdegree;
  stdoutdegree=sqrt(stdoutdegree);

  stdtotdegree/=(double)genecounter;
  stdtotdegree=stdtotdegree-avgtotdegree*avgtotdegree;
  stdtotdegree=sqrt(stdtotdegree);


  sprintf(fname,"%s/%s%s",despath,"DegreeDistr_",filename);
  f=fopen(fname,"w");
  for(i=0;i<15;i++)
  {
    fprintf(f,"%i\t%f\t%f\t%f\n",i,indegreedistr[i],outdegreedistr[i],totdegreedistr[i]);
    //fprintf(f,"%f\t%f\t%f\n",avgindegreedistr[i],avgoutdegreedistr[i],avgtotdegreedistr[i]);
  }
  fclose(f);

  /*for(i=0;i<15;i++)
  {
    avgindegreedistr[i]=0;
    avgoutdegreedistr[i]=0;
    avgtotdegreedistr[i]=0;
  }*/
  
}


void Agent::WriteGenomeAndNetworkProperties(int aid)
{
  FILE *f;
  char fname[800];

  sprintf(fname,"%s/%s",despath,"GenomeAndNetworkProperties");
  f=fopen(fname,"a");
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",G->glength_);//2
  fprintf(f,"%i\t",G->gnrgenes_);//3
  //fprintf(f,"%i\t",nridgenes);//4
  fprintf(f,"%i\t",G->gnrtfbs_);//5
  fprintf(f,"%i\t",nrplustfbs);//6
  fprintf(f,"%i\t",nrmintfbs);//7
  fprintf(f,"%f\t",avgindegree);//8
  fprintf(f,"%f\t",avgoutdegree);//9
  fprintf(f,"%f\t",avgtotdegree);//10
  fprintf(f,"%f\t",minindegree);//11
  fprintf(f,"%f\t",minoutdegree);//12
  fprintf(f,"%f\t",mintotdegree);//13
  fprintf(f,"%f\t",maxindegree);//14
  fprintf(f,"%f\t",maxoutdegree);//15
  fprintf(f,"%f\t",maxtotdegree);//16
  fprintf(f,"%f\t",stdindegree);//17
  fprintf(f,"%f\t",stdoutdegree);//18
  fprintf(f,"%f\t",stdtotdegree);//19
  fprintf(f,"\n");
  fclose(f);
}


void Agent::DetermineAttractorProperties()
{
  int t,i,j,k;
  int samefound;
  int nrtrials=100;//500000;//10000;
  int proteinlevels[5];
  double temp;
  static int agentcounter=0;
  FILE *f;
  char fname[800];

  //start with clean arrays
  nrequilibria=0;
  for(i=0;i<1000;i++)
    for(j=0;j<NrGeneTypes;j++)
      equilibria[i][j]=0;
  for(i=0;i<1000;i++)
    sizeequilibria[i]=0;

  //initialise proteinlevels array
  proteinlevels[0]=0;
  proteinlevels[1]=20;
  proteinlevels[2]=40;
  proteinlevels[3]=80;
  proteinlevels[4]=160;
  
  list<Cell>::iterator iter;
  iter=cells.begin();
  
  for(t=0;t<nrtrials;t++)
    {
      //zero startstate maternal genes
      for(i=0;i<NrMatGeneTypes;i++)
	(*iter).proteinstates[i]=0;
      //random startstate for proteinstates
      for(i=NrMatGeneTypes;i<NrGeneTypes;i++)
	{
	  //CL[0]->proteinstates[i]=uniform()*200;
	  temp=uniform();
	  if(temp<0.2)
	    (*iter).proteinstates[i]=proteinlevels[0];
	  else if(temp<0.4)
	    (*iter).proteinstates[i]=proteinlevels[1];
	  else if(temp<0.6)
	    (*iter).proteinstates[i]=proteinlevels[2];
	  else if(temp<0.8)
	    (*iter).proteinstates[i]=proteinlevels[3];
	  else
	    (*iter).proteinstates[i]=proteinlevels[4];
	}
      //startstate for signalstates
      for(i=0;i<NrSignGeneTypes/2;i++)//direct cell cell communication
	(*iter).signalstates[i]=0;//no input or yourself????
      for(i=NrSignGeneTypes/2;i<NrSignGeneTypes;i++)//diffusive cell cell communication
	(*iter).signalstates[i]=(*iter).proteinstates[i+NrMatGeneTypes];// your own conc as input

      //iteratively update cell state
      for(i=0;i<NrDevSteps;i++)
	{
	  N->UpdateNetworkState(1,(*iter).genestates,(*iter).proteinstates,(*iter).signalstates);

	  for(j=0;j<NrSignGeneTypes/2;j++)//direct cell cell communication
	    (*iter).signalstates[j]=0;//no input or yourself????

	  for(j=NrSignGeneTypes/2;j<NrSignGeneTypes;j++)//diffusive cell cell communication
	    (*iter).signalstates[j]=(*iter).proteinstates[j+NrMatGeneTypes];// your own conc as input
	}

      //determine whether new or old equilibrium
      samefound=0;
      for(i=0;i<nrequilibria;i++)//compare with other equilibria
	{
	  if(DifferentEquilibrium(equilibria[i],(*iter).proteinstates)==0)
	    {
	      samefound=1;
	      sizeequilibria[i]++;
	      break;
	    }
	}
      if(samefound==0)//new one
	{
	  for(j=0;j<NrGeneTypes;j++)
	    equilibria[nrequilibria][j]=(*iter).proteinstates[j];
	  nrequilibria++;
	}
    }


  //normalised sizes of attractors
  for(i=0;i<nrequilibria;i++)
    sizeequilibria[i]/=(double)nrtrials;
  //determine avg size of equilibria
  avgeqsize=0;
  for(i=0;i<nrequilibria;i++)
    avgeqsize+=sizeequilibria[i];
  if(nrequilibria>1)
    avgeqsize/=nrequilibria;
  
  //usage of equilibria during and at end of development
  int eqused[1000];
  int equsedatend[1000];
  nrequsedtotal=0;
  nrequsedend=0;
  for(i=0;i<1000;i++)
    {
      eqused[i]=0;
      equsedatend[i]=0;
    }
  for(i=0;i<NrStorages;i++)
    {
      for(j=0;j<NrFinalCells;j++)
	{
	  for(k=0;k<nrequilibria;k++)
	    {
	      if(DifferentEquilibrium((double *)E[i][j],equilibria[k])==0)
		{
		  eqused[k]++;
		  if(i==NrStorages-1)
		    equsedatend[k]++;
		  break;
		}
	    }
	}
    }
  for(i=0;i<nrequilibria;i++)
    {
      if(eqused[i]>0)
	nrequsedtotal++;
      if(equsedatend[i]>0)
	nrequsedend++;
    }

  if(agentcounter%100==0)
    {
      sprintf(fname,"%s/%s%i",despath,"Attractors",agentid);
      f=fopen(fname,"w");
      for(i=0;i<nrequilibria;i++)
	{
	  fprintf(f,"%i\t",i);
	  for(j=0;j<NrGeneTypes;j++)
	    fprintf(f,"%i\t",(int)equilibria[i][j]);
	  fprintf(f,"%f\t",sizeequilibria[i]);
	  fprintf(f,"\n");
	}
      fclose(f);
    }

  agentcounter++;
}



void Agent::WriteAttractorProperties(int aid)
{
  FILE *f;
  char fname[800];

  sprintf(fname,"%s/%s",despath,"AttractorProperties");
  f=fopen(fname,"a");
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",nrequilibria);//2
  fprintf(f,"%i\t",nrequsedtotal);//3
  fprintf(f,"%i\t",nrequsedend);//4
  fprintf(f,"%f\t",avgeqsize);//5
  fprintf(f,"%i\t",nrlonganddifbands);//6
  fprintf(f,"\n");
  fclose(f);
}


/*
 v oid Agent::OldDetermineL*oopAndMotifProperties()
{
  nrposloops=0;
  nrposloopsE=0;
  nrposloopsI=0;
  nrposmutual=0;
  nrposmutualE=0;
  nrposmutualI=0;
  nrnegmutual=0;
  nrnegmutualE=0;
  nrnegmutualI=0;
  nrposthreesome=0;
  nrposthreesomeE=0;
  nrposthreesomeI=0;

  //int A[200][200];
  int genetypes[200];
  int i,j;
  for(i=0;i<200;i++)
    for(j=0;j<200;j++)
      A[i][j]=0;
  for(i=0;i<200;i++)
    genetypes[i]=-1;
  
  
  int genecounter=0;  
  Genome::iter it;
  TFBS *tfbs;
  Gene *gene;

  //which genes are located where in genome
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  genetypes[genecounter]=gene->type;
	  genecounter++;
	}
      it++;
    }

  //which type of tfbs have genes upstream of them
  //multiple genes of the same type can bind to this tfbs
  int genecounter2=0;
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsTFBS(*it))//tfbs
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  for(i=0;i<genecounter;i++)
	    {
	      if(tfbs->type==genetypes[i])
		A[i][genecounter2]=tfbs->weight;
	    }
	}
      else
	genecounter2++;
      it++;
    }

  //determine positive autoloops
  for(i=0;i<genecounter;i++)
    {
      if(A[i][i]==1)
	{
	  nrposloops++;
	  if(i<NrMatGeneTypes+NrSignGeneTypes)
	    nrposloopsE++;
	  else 
	    nrposloopsI++;
	}
    }
  //determine plusplus and minmin mutual loops
  for(i=0;i<genecounter;i++)
   for(j=i+1;j<genecounter;j++)//to ensure non double counting
     {
       if(A[i][j]==1 && A[j][i]==1)
	 {
	   nrposmutual++;
	   if(i<NrMatGeneTypes+NrSignGeneTypes || j<NrMatGeneTypes+NrSignGeneTypes)
	     nrposmutualE++;
	   else 
	     nrposmutualI++;
	 }
       else if(A[i][j]==-1 && A[j][i]==-1)
	 {
	   nrnegmutual++;
	   if(i<NrMatGeneTypes+NrSignGeneTypes || j<NrMatGeneTypes+NrSignGeneTypes)
	     nrnegmutualE++;
	   else 
	     nrnegmutualI++;
	 }
     }

  
  int Adummy[200][200];
  for(i=0;i<200;i++)
    for(j=0;j<200;j++)
      Adummy[i][j]=A[i][j];

  int prod;
  int count;
  int I,J,K;
  int jj;
  for(i=0;i<genecounter;i++)
   for(j=i+1;j<genecounter;j++)//to ensure non double counting
     {
       if(Adummy[i][j]!=0)//potential start of a threesome loop
	 {
	   prod=Adummy[i][j];
	   count=1;
	   I=i;
	   J=j;
	   for(jj=0;jj<genecounter;jj++)
	     {
	       if(Adummy[J][jj]!=0)//potential second point on threesome loop
		 {
		   if(!(jj==i) && !(J==jj))//not a mutual loop or autoloop, so could still be threesome loop
		     {
		       prod*=Adummy[J][jj];
		       count=2;
		       K=jj;
		       if(Adummy[K][I]!=0)//third point on threesome loop
			 {
			   prod*=Adummy[K][I];
			   count=3;
			   
			   if(count==3 && prod==1)
			     {
			       nrposthreesome++;
			       if(I<NrMatGeneTypes+NrSignGeneTypes || 
				  J<NrMatGeneTypes+NrSignGeneTypes ||
				  K<NrMatGeneTypes+NrSignGeneTypes )
				 nrposthreesomeE++;
			       else
				 nrposthreesomeI++;
			     }
				 
			   
			   //restore previous setting and look for other second point
			   prod/=Adummy[K][I];
			   prod/=Adummy[J][jj];
			   count=1;
			 }
		       else
			 {
			   //restore previous setting and look for other second point
			   prod/=Adummy[J][jj];
			   count=1;
			 }
		     }
		 }
	     }
	 }
     }
}
  //which type of tfbs have genes upstream of them
  //multiple genes of the same type can bind to this tfbs
  int genecounter2=0;
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsTFBS(*it))//tfbs
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  for(i=0;i<genecounter;i++)
	    {
	      if(tfbs->type==genetypes[i])
		A[i][genecounter2]=tfbs->weight;
	    }
	}
      else
	genecounter2++;
      it++;
    }

  //determine positive autoloops
  for(i=0;i<genecounter;i++)
    {
      if(A[i][i]==1)
	{
	  nrposloops++;
	  if(i<NrMatGeneTypes+NrSignGeneTypes)
	    nrposloopsE++;
	  else 
	    nrposloopsI++;
	}
    }
  //determine plusplus and minmin mutual loops
  for(i=0;i<genecounter;i++)
   for(j=i+1;j<genecounter;j++)//to ensure non double counting
     {
       if(A[i][j]==1 && A[j][i]==1)
	 {
	   nrposmutual++;
	   if(i<NrMatGeneTypes+NrSignGeneTypes || j<NrMatGeneTypes+NrSignGeneTypes)
	     nrposmutualE++;
	   else 
	     nrposmutualI++;
	 }
       else if(A[i][j]==-1 && A[j][i]==-1)
	 {
	   nrnegmutual++;
	   if(i<NrMatGeneTypes+NrSignGeneTypes || j<NrMatGeneTypes+NrSignGeneTypes)
	     nrnegmutualE++;
	   else 
	     nrnegmutualI++;
	 }
     }

  
  int Adummy[200][200];
  for(i=0;i<200;i++)
    for(j=0;j<200;j++)
      Adummy[i][j]=A[i][j];

  int prod;
  int count;
  int I,J,K;
  int jj;
  for(i=0;i<genecounter;i++)
   for(j=i+1;j<genecounter;j++)//to ensure non double counting
     {
       if(Adummy[i][j]!=0)//potential start of a threesome loop
	 {
	   prod=Adummy[i][j];
	   count=1;
	   I=i;
	   J=j;
	   for(jj=0;jj<genecounter;jj++)
	     {
	       if(Adummy[J][jj]!=0)//potential second point on threesome loop
		 {
		   if(!(jj==i) && !(J==jj))//not a mutual loop or autoloop, so could still be threesome loop
		     {
		       prod*=Adummy[J][jj];
		       count=2;
		       K=jj;
		       if(Adummy[K][I]!=0)//third point on threesome loop
			 {
			   prod*=Adummy[K][I];
			   count=3;
			   
			   if(count==3 && prod==1)
			     {
			       nrposthreesome++;
			       if(I<NrMatGeneTypes+NrSignGeneTypes || 
				  J<NrMatGeneTypes+NrSignGeneTypes ||
				  K<NrMatGeneTypes+NrSignGeneTypes )
				 nrposthreesomeE++;
			       else
				 nrposthreesomeI++;
			     }
				 
			   
			   //restore previous setting and look for other second point
			   prod/=Adummy[K][I];
			   prod/=Adummy[J][jj];
			   count=1;
			 }
		       else
			 {
			   //restore previous setting and look for other second point
			   prod/=Adummy[J][jj];
			   count=1;
			 }
		     }
		 }
	     }
	 }
     }
}*/

bool compare_ints (int first,int second)
{
  return ( first < second);
}


/*
void Agent::removeDoubles(list <int> &tallying, int finalcall)
{
   
  list< list <int> > counterstore; //required only to find the double-counted loops
  list <int> :: iterator il, el;
  list<int>::iterator tal, ntal;

  list< list<int> >::iterator walk, next, cwalk, cnext;
  
  int ns=0;
  list< list<int> > copy;
  
  copy=loopstore;
  printf("size of loopstore: %d\n", (int)loopstore.size());
  
  for(walk=loopstore.begin(), tal=tallying.begin(), cwalk=copy.begin(); walk!=loopstore.end(); ++walk, ++tal, ++cwalk)
  {
    next=walk;
    next++;
    cnext=cwalk;
    cnext++;
    ntal=tal;
    ntal++;
    
    while(next!=loopstore.end()){
      ns=0;
      
      if ((int)(*walk).size()==(int)(*next).size())
      {
	//printf("comparing...\n");
	(*cwalk).sort(compare_ints);
	(*cnext).sort(compare_ints);
	for(il=(*cwalk).begin(),el=(*cnext).begin(); il!=(*cwalk).end(); ++il, ++el)
	{
	  //printf("sorted: %d\n",(*il));
	  if((*il)!=(*el))
	  {
	    ns=1;
	    break;
	  }
	}
	if(!ns){
	  next=loopstore.erase(next);
	  cnext=copy.erase(cnext);
	  ntal=tallying.erase(ntal);
	  continue;
	}
      }
      ++next;
      ++cnext;
      ++ntal;
    }
    if (finalcall){
      if((*tal)>0)
	nrpos++;
      else
	nrneg++;
    }
  }
}
void Agent::DetermineLoopAndMotifProperties()
{
  nrposauto=0;
  nrnegauto=0;
  nrpos=0;
  nrneg=0;
  int A[200][200];

  compreg=0; //number of complex (+-) regulations
  
  //int A[200][200];
  int genetypes[200];
  int i,j;
  for(i=0;i<200;i++)
    for(j=0;j<200;j++)
      A[i][j]=0;
  for(i=0;i<200;i++)
    genetypes[i]=-1;
  
  
  int genecounter=0;  
  Genome::iter it;
  TFBS *tfbs;
  Gene *gene;

  //which genes are located where in genome
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  genetypes[genecounter]=gene->type;
	  genecounter++;
	}
      it++;
    }

  //which type of tfbs have genes upstream of them
  //multiple genes of the same type can bind to this tfbs
  int genecounter2=0;
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
  {
    if(G->IsTFBS(*it))//tfbs
    {
      tfbs=dynamic_cast<TFBS *>(*it);
      for(i=0;i<genecounter;i++)
      {
    
	if(tfbs->type==genetypes[i])
	{
	  if(genetypes[i]==genetypes[genecounter2]) //autoloop: don't put in matrix
            if(tfbs->weight<0)
	      nrnegauto++;
	    else
	      nrposauto++;

	  else if((tfbs->weight<0)==(A[i][genecounter2]<0) || !A[i][genecounter2]) //same sign or A is still empty
		A[i][genecounter2]+=tfbs->weight;
	  
	  else
	  {
	    A[i][genecounter2]=1; //a choice has to be made here: A is integer, so cannot account for having both neg and pos interaction. we'll count it as an "on" interaction for now.
	    compreg++;
	  }
	}
      }
    }
    else
      genecounter2++;
    it++;
  }
  
 
  
  
  // int Adummy[200][200];
//  for(i=0;i<200;i++)
//    { 
//    for(j=0;j<200;j++)
//      if(A[i][j])
//        printf("gene %d %d gene %d\n",genetypes[i], A[i][j], genetypes[j]);
// }
  
  list <int> nodes, counter, tallying; //temporary storage for the genetypes in the path followed, their genecounter (position in genome, unique identifier), and to keep score whether neg or pos loop.
  int genenr; 
  int tally=1; //to find the 
 
  list< list <int> > counterstore; //required only to find the double-counted loops
  list <int> :: iterator il, el;
   
  for(int row=0; row<genecounter; row++)
  {
    //empty the temporary storage
    nodes.clear();
    counter.clear();
    //push first gene in list
    nodes.push_back(genetypes[row]);
    counter.push_back(row);
    tally=1;
    genenr=row;
    printf("start gene %d, ",nodes.back());
    printf("row = %d\n", row);
    
    for(int col=0; col<genecounter; col++) //go through all genes that gene "genenr" regulates
    {
      if(A[genenr][col] && find(counter.begin(),counter.end(), col)==counter.end()){ //gene "genenr" regulates this gene, we have no loop yet
	//push gene at end of list
	nodes.push_back(genetypes[col]);
	counter.push_back(col);
	tally*=A[genenr][col];
	genenr=col;
	col=-1;
	//printf("tally: %d\n",tally);
	//printf("next gene %d, counter=%d ", nodes.back(), counter.back());
      }
      else if (A[genenr][col] && find(counter.begin(), counter.end(), col)!=counter.end()) //found a loooop
      {
	il=find(nodes.begin(), nodes.end(), genetypes[col]); //find same position in storage of the genes (flawed: may be multipe genes of same type in a loop)
	//printf("loop to %d\n",(*il));
	list <int> vec (il, nodes.end());
	tallying.push_back(tally);
	//printf("el %d\n",(*vec.begin()));
	loopstore.push_back(vec);
	genenr=counter.back();
	//printf("genenr= %d\n", genenr);
      }
      
      while(col==genecounter-1 && nodes.size()>1) 
      {
	nodes.pop_back();
	//genenr=nodes.back();
	col=counter.back();
	counter.pop_back();
	genenr=counter.back();
	//printf("back to gene %d, counter %d\n", nodes.back(), col);
      }
    }
    
    //if (!(row%2)){
      printf("am here!!\n");
      removeDoubles(tallying,0);
    //}
  }
  
  removeDoubles(tallying,1);
  
  
  
  
}*/

void Agent::removeDoubles_type()
{
  nrpos=0;
  nrneg=0;
  
  vector<int> :: iterator il, el;
  vector<int>::iterator tal, ntal;

  list< vector<int> >::iterator walk, next, cwalk, cnext, pwalk, pnext; //loopstore, copy of loopstore, counterstore
  
  int ns=0;
  list< vector<int> > copy;
  
  copy=loopstore;
   
  for(walk=loopstore.begin(), tal=tallying.begin(), cwalk=copy.begin(), pwalk=counterstore.begin(); walk!=loopstore.end(); ++walk, ++tal, ++cwalk, ++pwalk)
  {
    next=walk;
    next++;
    cnext=cwalk;
    cnext++;
    pnext=pwalk;
    pnext++;
    ntal=tal;
    ntal++;
    
    while(next!=loopstore.end()){
      ns=0;
      
      if ((int)(*walk).size()==(int)(*next).size() && (*tal)==(*ntal) )
      {
	sort((*cwalk).begin(),(*cwalk).end());
	sort((*cnext).begin(),(*cnext).end());
	for(il=(*cwalk).begin(),el=(*cnext).begin(); il!=(*cwalk).end(); ++il, ++el)
	{
	  //printf("sorted: %d\n",(*il));
	  if((*il)!=(*el))
	  {
	    ns=1;
	    break;
	  }
	}
	if(!ns){ //the two vectors are the same
	  next=loopstore.erase(next);
	  cnext=copy.erase(cnext);
	  pnext=counterstore.erase(pnext);
	  ntal=tallying.erase(ntal);
	  continue; //don't have to move the iterator anymore
	}
      }
      ++next;
      ++cnext;
      ++ntal;
      ++pnext;
    }
    
    if((*tal)>0)
      nrpos++;
    else
      nrneg++;
  }//loop through loopstore
  
}

void Agent::removeDoubles_id(int finalcall)
{
//  removing double counted loops from loopstore, tallying the number of positive and negative loops 
  
  vector<int>::iterator tal, ntal, il, el;
  list< vector<int> >::iterator walk, next, cwalk, cnext;
  
  int ns=0;
  //list< vector<int> > copy;
  
  //copy=loopstore;
 
  for(walk=loopstore.begin(), tal=tallying.begin(), cwalk=counterstore.begin(); walk!=loopstore.end(); ++walk, ++tal, ++cwalk)
  {
    next=walk;
    next++;
    cnext=cwalk;
    cnext++;
    ntal=tal;
    ntal++;
    
    while(next!=loopstore.end()){
      ns=0;
      
      if ((int)(*walk).size()==(int)(*next).size() && (*tal)==(*ntal)) //same length and same sign
      {
	//printf("comparing...\n");
	sort((*cwalk).begin(),(*cwalk).end());
	sort((*cnext).begin(),(*cnext).end());
	
	for(il=(*cwalk).begin(),el=(*cnext).begin(); il!=(*cwalk).end(); ++il, ++el)
	{
	  //printf("sorted: %d\n",(*il));
	  if((*il)!=(*el))
	  {
	    ns=1;
	    break;
	  }
	}
	if(!ns){
	  next=loopstore.erase(next);
	  cnext=counterstore.erase(cnext);
	  ntal=tallying.erase(ntal);
	  continue;
	}
      }
      ++next;
      ++cnext;
      ++ntal;
    }
    if (finalcall){
      if((*tal)>0)
	nrpos++;
      else
      nrneg++;
    }
  }
}

void Agent::DetermineLoopAndMotifProperties()
{
  nrposauto=0;
  nrnegauto=0;
  nrpos=0;
  nrneg=0;
  int A[200][200];

  compreg=0; //number of complex (+-) regulations
  //printf("am here!\n");
  //int A[200][200];
  int genetypes[200];
  int i,j;
  for(i=0;i<200;i++)
    for(j=0;j<200;j++)
      A[i][j]=0;
  for(i=0;i<200;i++)
    genetypes[i]=-1;
  
  
  int genecounter=0;  
  Genome::iter it;
  TFBS *tfbs;
  Gene *gene;

  //which genes are located where in genome
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  genetypes[genecounter]=gene->type;
	  genecounter++;
	}
      it++;
    }

  //which type of tfbs have genes upstream of them
  //multiple genes of the same type can bind to this tfbs
  int genecounter2=0;
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
  {
    if(G->IsTFBS(*it))//tfbs
    {
      tfbs=dynamic_cast<TFBS *>(*it);
      for(i=0;i<genecounter;i++)
      {
    
	if(tfbs->type==genetypes[i])
	{
	  if(genetypes[i]==genetypes[genecounter2]) //autoloop: don't put in matrix
            if(tfbs->weight<0)
	      nrnegauto++;
	    else
	      nrposauto++;

	  else if((tfbs->weight<0)==(A[i][genecounter2]<0) || !A[i][genecounter2]) //same sign or A is still empty
		A[i][genecounter2]+=tfbs->weight;
	  
	  else
	  {
	    A[i][genecounter2]=1; //a choice has to be made here: A is integer, so cannot account for having both neg and pos interaction. we'll count it as an "on" interaction for now.
	    compreg++;
	  }
	}
      }
    }
    else
      genecounter2++;
    it++;
  }
  
 
  
  //debugging
 /* for(i=0;i<200;i++)
  { 
    for(j=0;j<200;j++)
      if(A[i][j])
	printf("gene %d %d gene %d\n",genetypes[i], A[i][j], genetypes[j]);
  }*/
  
//    go through interaction matrix, find the loops 
  
  vector <int> nodes, counter; //temporary storage for the genetypes in the path followed, their genecounter (position in genome, unique identifier), and to keep score whether neg or pos loop.
  int genenr; 
  int tally=1; //to find whether loop is positive or negative
  
  int index=0; //to find position of iterator
  int found=0;
  
  int loopsize=0; //keep track of loop size: break free if it becomes too big... NEW!!
  
  vector <int> :: iterator il, el;
   
  for(int row=0; row<genecounter; row++)
  {
    //empty the temporary storage
    nodes.clear();
    counter.clear();
    //push first gene in list
    nodes.push_back(genetypes[row]);
    counter.push_back(row);
    tally=1;
    genenr=row;
    loopsize=1;
    //printf("start gene %d, ",nodes.back());
    //printf("row = %d\n", row);
    
    for(int col=0; col<genecounter; col++) //go through all genes that gene "genenr" regulates
    {
      if(A[genenr][col] && find(counter.begin(),counter.end(), col)==counter.end()){ //gene "genenr" regulates this gene, we have no loop yet
	//push gene at end of list
	nodes.push_back(genetypes[col]);
	counter.push_back(col);
	//tally*=A[genenr][col];
	genenr=col;
	col=-1;
	loopsize++;
	//printf("tally: %d\n",tally);
	//printf("next gene %d, counter=%d ", nodes.back(), counter.back());
      }
      else if (A[genenr][col] && (el=find(counter.begin(), counter.end(), col))!=counter.end()) //found a loooop
      {
	//to find the same position in the nodes vector
	index=el-counter.begin();
	il=nodes.begin()+index;
	//tally*=A[genenr][col];
	vector <int> vec (il, nodes.end());
	vector <int> cvec (el, counter.end());
	//store whether pos or neg fb loop, and the loop itself both as gene type and the unique identifier.
	//tallying.push_back(tally);
	counterstore.push_back(cvec);
	loopstore.push_back(vec);
	genenr=counter.back();
	found=1;
	
	for (int i=0; i<vec.size(); i++)
	{
	  //printf("gene %d, id %d\t",vec[i],cvec[i]);
	  if(i<vec.size()-1) tally*=A[cvec[i]][cvec[i+1]];
	  else tally*=A[cvec[i]][cvec[0]]; 
	  
	}

	tallying.push_back(tally);
	tally=1;
      }
      
            
      while( (found || loopsize >15 || (!A[row][col]  && col==genecounter-1) ) && nodes.size()>1 ) //iterate back to see if we can find more loops
      {
	//if(A[row][col]){
	  //printf("changing tally; was %d, is now %d; row=%d, col=%d\n",tally,tally/A[row][col], row, col );
	  //tally=tally/A[row][col];
	//}
	nodes.pop_back();
	//genenr=nodes.back();
	col=counter.back();
	counter.pop_back();
	genenr=counter.back();
	loopsize--;
	found=0; 
	
	//printf("back to gene %d, counter %d\n", nodes.back(), col);
      }
      
    }//col loop
    
    removeDoubles_id(0);
  } //row loop 
 
  removeDoubles_id(1);
}

/*void Agent::WriteLoopAndMotifProperties(int aid)
{
  FILE *f;
  char fname[800];

  sprintf(fname,"%s/%s",despath,"LoopAndMotifProperties");
  f=fopen(fname,"a");
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",nrposloops);//2
  fprintf(f,"%i\t",nrposmutual);//5
  fprintf(f,"%i\t",nrnegmutual);//8
  fprintf(f,"%i\t",nrnegmutualI);//9
  fprintf(f,"%i\t",nrnegmutualE);//10
  fprintf(f,"%i\t",nrposthreesome);//11
  fprintf(f,"%i\t",nrposthreesomeI);//12
  fprintf(f,"%i\t",nrposthreesomeE);//13
  fprintf(f,"%i\t",nrposloops+nrposmutual+nrnegmutual+nrposthreesome);//14
  fprintf(f,"%i\t",nrposloopsI+nrposmutualI+nrnegmutualI+nrposthreesomeI);//15
  fprintf(f,"%i\t",nrposloopsE+nrposmutualE+nrnegmutualE+nrposthreesomeE);//16
  fprintf(f,"%i\t",nrposloops+nrposmutual+nrnegmutual);//17
  fprintf(f,"%i\t",nrposloopsI+nrposmutualI+nrnegmutualI);//18
  fprintf(f,"%i\t",nrposloopsE+nrposmutualE+nrnegmutualE);//19
  fprintf(f,"\n");
  fclose(f);
}*/

void Agent::WriteLoopAndMotifProperties(int aid, char *name)
{
  FILE *f, *f2;
  char fname[800];
  int countarray[20];
    
  for(int i=0;i<20;i++)
    countarray[i]=0;
    
  sprintf(fname,"%s/%s_%s",despath,"LoopAndMotifProperties", name);
  f=fopen(fname,"w");
  
 list<vector<int> >::iterator walk; 
 vector<int>::iterator el; 
  

  /** some generic properties **/
  
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",nrposauto);//2
  fprintf(f,"%i\t",nrnegauto);//3
  fprintf(f,"%d\t", nrpos); //4
  fprintf(f,"%d\t", nrneg);//5
  removeDoubles_type();
  fprintf(f,"%d\t", nrpos); //6 //the non-redundant loops
  fprintf(f,"%d\n", nrneg);//7
  fprintf(f, "\n\n");
  //fprintf(f,"\n\nLoops of size:");
  /** print the loops themselves to a separate file **/
  for(walk=loopstore.begin(); walk!=loopstore.end(); ++walk)
  {
    if((*walk).size()<=15)
      countarray[(*walk).size()]++; //tally loop sizes
      else
	countarray[16]++;
      
    for(el=(*walk).begin(); el!=(*walk).end();++el) //print loop
   {
     fprintf(f,"%d\t",(*el));
   }
   fprintf(f,"\n\n");
  }
  
  fclose(f);
  
  /** ****** */
  sprintf(fname,"%s/%s_%s",despath,"LoopSizeFreqs", name); //corrected for doubles in gene type! 
  f2=fopen(fname,"w");
  
  for(int i=0; i<16;i++)
    fprintf(f2,"%d\t%d\n", i, countarray[i]); 
  
  fclose(f2);
  
}


void Agent::DetermineGenoNetworkPhenoDistance(Agent *M,double *genodist,double *phenodist)
{
  int i,j;
  double gdist=0;
  double pdist=0;
  int geneocc1[NrGeneTypes];
  int geneocc2[NrGeneTypes];
  int tfbsocc1[NrGeneTypes][NrGeneTypes*2];
  int tfbsocc2[NrGeneTypes][NrGeneTypes*2];
  Network::iterv iv;
  Network::iterel el;
  Network::itere e;
  int gene_genetype;
  int tfbs_genetype;
  int weight;

  //genotypic distance
  //init
  for(i=0;i<NrGeneTypes;i++)
    {
      geneocc1[i]=0;
      geneocc2[i]=0;
    }
  for(i=0;i<NrGeneTypes;i++)
    for(j=0;j<2*NrGeneTypes;j++)
      {
	tfbsocc1[i][j]=0;
	tfbsocc2[i][j]=0;
      }
  //genome 1
  iv=N->VL->begin();
  el=N->AL->begin();
  while(iv!=N->VL->end())
    {
      gene_genetype=(*iv)->Gen->type;
      geneocc1[gene_genetype]++;
      for(e=el->begin();e!=el->end();e++)
	{
	  tfbs_genetype=(*e)->V->Gen->type;
	  weight=(*e)->weight;
	  if(weight==-1)
	    tfbsocc1[gene_genetype][tfbs_genetype]++;
	  else
	    tfbsocc1[gene_genetype][12+tfbs_genetype]++;	      
	}
      iv++;
      el++;
    }
  //genome 2
  iv=M->N->VL->begin();
  el=M->N->AL->begin();
  while(iv!=M->N->VL->end())
    {
      gene_genetype=(*iv)->Gen->type;
      geneocc2[gene_genetype]++;
      for(e=el->begin();e!=el->end();e++)
	{
	  tfbs_genetype=(*e)->V->Gen->type;
	  weight=(*e)->weight;
	  if(weight==-1)
	    tfbsocc2[gene_genetype][tfbs_genetype]++;
	  else
	    tfbsocc2[gene_genetype][12+tfbs_genetype]++;	      
	}
      iv++;
      el++;
    }
  //distance
  for(i=0;i<NrGeneTypes;i++)
    gdist+=abs(geneocc1[i]-geneocc2[i]);
  for(i=0;i<NrGeneTypes;i++)
    for(j=0;j<2*NrGeneTypes;j++)
      gdist+=abs(tfbsocc1[i][j]-tfbsocc2[i][j]);
  gdist/=(double)(NrGeneTypes);

  (*genodist)=gdist;


  //phenotypic distance: compare celltypes
  pdist=0;
  for(i=0;i<NrFinalCells;i++)
    {
      if(types[NrDevSteps-1][i]!=M->types[NrDevSteps-1][i])
	pdist++;
    }
  pdist/=(double)NrFinalCells;

  (*phenodist)=pdist;
}

void Agent::WriteGenoNetworkPhenoDistance(int aid,double genodist,double phenodist,double cumgenodist,double cumphenodist)
{
  //indien zelfde als vorige agent: afstanden nul, dus cumulatieve afstanden constant
  
  FILE *f;
  char fname[800];

  sprintf(fname,"%s/%s",despath,"GenoPhenoDistances");
  f=fopen(fname,"a");
  fprintf(f,"%i\t",agentid);//1
  fprintf(f,"%f\t",genodist);//2
  fprintf(f,"%f\t",phenodist);//3
  fprintf(f,"%f\t",cumgenodist);//4
  fprintf(f,"%f\t",cumphenodist);//5
  fprintf(f,"%i\t",nrlonganddifbands);//6
  fprintf(f,"\n");
  fclose(f);
}