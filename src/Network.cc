#include "Network.hh"


Network::Network()
{
  VL=NULL;
  AL=NULL;
}

Network::~Network()
{  
  iterv iv;
  iterel el;
  itere e;

  if(VL!=NULL)
    {
      //delete what iterators point to
      for(iv=VL->begin();iv!=VL->end();iv++)
	{
	  delete (*iv);
	}

      //delete list of pointers itself
      iv=VL->erase(VL->begin(),VL->end());
      //delete pointer to list
      delete VL;
      VL=NULL;
    }

  if(AL!=NULL)
    {
      for(el=AL->begin();el!=AL->end();el++)
	{
	  for(e=el->begin();e!=el->end();e++)
	    {
	      delete (*e);
	    }
	  e=el->erase(el->begin(),el->end());
	}
      //delete list of pointers itself
      el=AL->erase(AL->begin(),AL->end());
      //delete pointer to list
      delete AL;
      AL=NULL;
    }
}


void Network::BuildNetwork(Genome *G)
{
  Genome::iter i;
  Genome::reviter ri;
  iterv iv;
  iterel el;
  itere e;
  int genecounter;
  Gene *gene;
  TFBS *tfbs;
  Vertex *ve;
  Edge *ed;
  int typematch;

  //first pass through genome:
  //build list of vertices
  VL=new list<Vertex *>();
  genecounter=0;
  for(i=(*G).ChromBBList->begin();i!=(*G).ChromBBList->end();i++)//forward
    {
      if((*G).IsGene((*i)))
	{
	  gene=dynamic_cast<Gene *>(*i);
	  ve=new Vertex(gene,genecounter);
	  (*VL).push_back(ve);
	  genecounter++;
	}
    }
 
  int first;
 
  //second pass through genome
  //build adjacency list
  AL=new list<list <Edge *> >();
  for(ri=(*G).ChromBBList->rbegin();ri!=(*G).ChromBBList->rend();ri++)//backward!!
    {
      if((*G).IsGene((*ri)))
	 {
	   AL->push_front( list<Edge *>() );//start new sublist for new gene
	   gene=dynamic_cast<Gene *>(*ri);
	 }
       else
	 {
	   tfbs=dynamic_cast<TFBS *>(*ri);//get TFBS
	   typematch=(*tfbs).type;//get to which genes it matches, may be more than 1
	   
	   first=0;
	   for(iv=(*VL).begin();iv!=(*VL).end();iv++)//look in VL not genome!
	     {
	       ve=(*iv);//get vertex
	       gene=((*ve).Gen);//find gene from vertex
	       if(typematch==(*gene).type)//see if gene matches tfbsite
		 {	   
		   if(first==0)//first gene corresponding to that binding site
		     {
		       ed=new Edge(ve,(*tfbs).weight,1,(*tfbs).HH);
		       AL->front().push_front(ed);
		       first=1;
		     }
		   else//more genes corresponding to that binding site
		     {
		       ed=new Edge(ve,(*tfbs).weight,0,(*tfbs).HH);
		       AL->front().push_front(ed);
		     }
		 }
	     }
	 }
    }
}


void Network::UpdateNetworkState(int steps,double genestates[300], double proteinstates[NrGeneTypes],double signalstates[NrSignGeneTypes])
{
  int i,j;
  iterv iv;
  iterel el;
  itere e;
  double transcract;
  double maxtranscract;
  double transcrrepr;
  double trscr;
  double efftranscr;
  double efftranscr2;
  double enhancer;
  int protein;
  double proteinstate;
  double Hstate;
  double transcr[NrGeneTypes];
  int genetype;
  double mu;
  double sigma;
  double noiseenhancer;
  int genecounter;

  for(i=0;i<steps;i++)
    {
      for(j=0;j<NrGeneTypes;j++)
	transcr[j]=0;

      iv=(*VL).begin();
      el=(*AL).begin();
      
      genecounter=0;
      //first step:
      //update gene expression state based on current gene states 
      while(iv!=(*VL).end())
	{
	  #ifndef FREEMORPH
	  if((*iv)->Gen->type<NrMatGeneTypes)
	    ;//not influenced, externally determined
	  else
	  {
	    #endif
	    #ifdef SUMINTEGRATION
	    trscr=0;
	    #endif
	    #ifdef MULTIPLYINTEGRATION
	    transcract=0;
	    maxtranscract=0;
	    transcrrepr=1;  
	    #endif
	    
	    for(e=el->begin();e!=el->end();e++)
	    {
	      if((*e)->use==1)
	      {            
		//use only first gene corresponding to certain TFBS
		//genes of same type are already added up in protein conc
		//otherwise you would use them double, added up in protein conc
		//and link by link
		protein=(*e)->V->Gen->type;
		if((*e)->V->Gen->type < NrMatGeneTypes)//maternal TF
		  proteinstate=__gnu_cxx::power(proteinstates[protein],N);
		else if((*e)->V->Gen->type < NrMatGeneTypes+NrSignGeneTypes)
		  proteinstate=__gnu_cxx::power(signalstates[protein-NrMatGeneTypes],N);//cell cel TF
		else
		  proteinstate=__gnu_cxx::power(proteinstates[protein],N);//normal TF
		Hstate=__gnu_cxx::power((*e)->HH,N);
		
		#ifdef SUMINTEGRATION
		trscr+=(*e)->weight*(proteinstate/(Hstate+proteinstate));
		#endif
		#ifdef MULTIPLYINTEGRATION
		if((*e)->weight>0)//activate gene expression
		{
		  transcract=proteinstate/(Hstate+proteinstate);
		  if(transcract>maxtranscract)
		    maxtranscract=transcract;
		}
		else if((*e)->weight<0)//repress gene expression
		  transcrrepr*=Hstate/(Hstate+proteinstate);
		#endif
		  		
	      }
	    }//end for loop over edges coming in on this vertex of network
	    //transcription of this gene contributes to expression of 
	    //protein type it codes for
	    
	    #ifdef SUMINTEGRATION
	    efftranscr=trscr;//-(*iv)->Gen->Th;
	    if(efftranscr>0)
	    {
	      efftranscr2=efftranscr*efftranscr;
	      enhancer=(efftranscr2/(efftranscr2+1))*(*iv)->Gen->EE;
	    }
	    else
	      enhancer=0;
	    #endif
	    
	    #ifdef MULTIPLYINTEGRATION
	    enhancer=maxtranscract*transcrrepr*(*iv)->Gen->EE;
	    #endif
	    
	    genetype=(*iv)->Gen->type;
	    transcr[genetype]+=enhancer;

	    //genetype=(*iv)->Gen->type;
	    genestates[genecounter]+=RungeKutta4(enhancer, (*iv)->Gen->DD,genestates[genecounter]);
	    //transcr[genetype]+=enhancer;

	    #ifndef FREEMORPH  
	  }
	  #endif
	  //after being finished with all edges of this one gene/vertex, go to next:
	  iv++;
	  el++;	  
	  genecounter++;
	}//end while loop over all vertices (genes) network    
      
      genecounter=0;
      //updating of gene states
      for(i=0;i<NrGeneTypes;i++)
      {
	#ifndef FREEMORPH
	if(i<NrMatGeneTypes)
	  continue;
	#endif
	proteinstates[i]=0;  
      }
      iv=(*VL).begin();
      while(iv!=(*VL).end())
      {
	genetype=(*iv)->Gen->type;
	
	#ifndef FREEMORPH
	if(genetype>=NrMatGeneTypes)
	{
	  #endif
	  proteinstates[genetype]+=genestates[genecounter];
	  
	  if((int)(proteinstates[i])<0 ||(int)(proteinstates[i])>25000)
	  {
	    printf("protein conc outside range %f\n",proteinstates[i]);
	    exit(1);
	  }
	  #ifndef FREEMORPH
	}
	#endif
	genecounter++;
	iv++;
      }
    }//end for time iteration loop
}//end function


