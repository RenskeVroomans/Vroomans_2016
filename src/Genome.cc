#include "Genome.hh"
#include "Misc.hh"

Genome::Genome()
{
  int i;
  ChromBBList=NULL;
}

Genome::~Genome()
{
  iter i;

  if(ChromBBList!=NULL)
    {  
      i=ChromBBList->begin();
      while(i!=ChromBBList->end())
	{
	  delete (*i);
	  i++;
	}
      i=ChromBBList->erase(ChromBBList->begin(),ChromBBList->end());
      delete ChromBBList;
      ChromBBList=NULL;
    }
}


void Genome::CopyPartOfGenome(iter begin,iter end)
{
  iter ii;
  ChromBB *cbb;

  ii=begin;
  while(ii!=end)
    {
      cbb=(*ii)->clone();
      (*ChromBBList).push_back(cbb);
      ii++;
    }
}


void Genome::CopyPartOfGenomeToTempList(iter begin,iter end,list<ChromBB*> &ChromBBListTemp)
{
  iter ii;
  ChromBB *cbb;
  
  ii=begin;
  while(ii!=end)
    {
      cbb=(*ii)->clone();
      ChromBBListTemp.push_back(cbb);
      ii++;
    }
}



void Genome::CloneGenome(const Genome *c)
{
  int i;

  ChromBBList=new list<ChromBB*>();

  CopyPartOfGenome(c->ChromBBList->begin(),c->ChromBBList->end());
  glength_=c->glength_;
  gnrgenes_=c->gnrgenes_;
  gnrtfbs_=c->gnrtfbs_;
  
}

void Genome::GenerateGenome()
{
  int i,j;
  Gene *gene;
  TFBS *tfbs;
  static int genesused[NrGeneTypes];
  int genenr;
  int ran;
  double temp;

  ChromBBList=new list<ChromBB*>();

  int nrcon;
  int type;

  glength_=0;
  gnrgenes_=0;
  gnrtfbs_=0;

  for(i=0;i<NrGeneTypes;i++)
    genesused[i]=0;

  for(i=0;i<NrGeneTypes;i++)
    {
      
      genenr=-1;
      while(genenr==-1)
	{
	  ran=(int)(uniform()*NrGeneTypes);
	  if(genesused[ran]==0)
	    {
	      genenr=ran;
	      genesused[ran]=1;
	    }
	}

      if(genenr<NrMatGeneTypes)
	;//no upstream region, conc externally determined
      else
	{
	  nrcon=1+(int)(uniform()*(AvCon+1));
	  for(j=0;j<nrcon;j++)
	    {
	      type=(int)(uniform()*NrGeneTypes);
	      
	      if(uniform()<0.5)
		tfbs=new TFBS(type,+1);
	      else
		tfbs=new TFBS(type,-1);
	      
	      glength_++;
	      gnrtfbs_++;
	      (*ChromBBList).push_back(tfbs);
	    }
	}
      gene=new Gene(genenr);
      glength_++;
      gnrgenes_++;
      (*ChromBBList).push_back(gene);
    }
  ListContent();
}


void Genome::CreateGenomeFromFile(char *fname, int agentidnr)
{
  FILE *f;

  //int i,j;
  Gene *gene;
  TFBS *tfbs;
  //static int genesused[NrGeneTypes];
  int type;
  int elementtype;
  //double Emax;
  int genetypecounter;
  char thisfile[1000];
  
  ChromBBList=new list<ChromBB*>();
  glength_=0;
  gnrgenes_=0;
  gnrtfbs_=0;
  
  if(!agentidnr){
    f=fopen(fname,"r");
    genetypecounter=0;
  }
  else{
    sprintf(thisfile,"%s/CodedGenomeAgent%.10d",fname,agentidnr);
    f=fopen(thisfile,"r");
    genetypecounter=0;
  }
  
  while(fscanf(f,"%i",&elementtype)!=EOF)
    {
      if(elementtype==-1)//binding site for repressing tf
	{
	  fscanf(f,"%i",&type);
	  tfbs=new TFBS(type,elementtype);
	  glength_++;
	  gnrtfbs_++;
	  (*ChromBBList).push_back(tfbs);
	  
	}
      else if(elementtype==1)//binding site for activating tf
	{
	  fscanf(f,"%i",&type);
	  tfbs=new TFBS(type,elementtype);
	  glength_++;
	  gnrtfbs_++;
	  (*ChromBBList).push_back(tfbs);
	}
      else if(elementtype==0)//gene
	{
	  fscanf(f,"%i",&type);
	  gene=new Gene(type);
	  glength_++;
	  gnrgenes_++;
	  (*ChromBBList).push_back(gene);
	}
    }
  fclose(f);

  //printf("created genome\n");
}

bool Genome::IsGene(ChromBB *cbb) const
{
  return (bool)(typeid(*cbb) == typeid(Gene));
}

bool Genome::IsTFBS(ChromBB *cbb) const
{
  return (bool)(typeid(*cbb) == typeid(TFBS));
}


int Genome::MutateGenome(int time)
{
  iter i;
 
  i=(*ChromBBList).begin();
  while(i!=(*ChromBBList).end())
    {
      if(IsGene(*i))
	{
	  i=GeneMutate(i);
	}
      else if(IsTFBS(*i))
	{
	  i=TFBSMutate(i);
	}
      else
	{
	  i++;
	}
     }
  TFBSInnovation(); 

  return 0;
}


void Genome::TFBSInnovation()
{
  int type;
  iter jj;
  TFBS *tfbsnew;

  if(uniform()<probtfbsinnov*glength_)
    {
      type=(int)(uniform()*NrGeneTypes);
      if(uniform()<0.5)
	tfbsnew=new TFBS(type,1);
      else
	tfbsnew=new TFBS(type,-1);
      jj=FindRandomGenePosition();
      jj=(*ChromBBList).insert(jj,tfbsnew);
      gnrtfbs_++;
      glength_++;
    }
}

Genome::iter Genome::FindFirstTFBSInFrontOfGene(iter ii) const
{
  reviter rii(ii);//init rev iter with passed iter
                  //note that it automatically
                  //points to the element BEHIND the
                  //element pointed to by iter
  reviter jj=(*ChromBBList).rend();//search should be bounded
                                   //by beginning of genome
                                   //which is rend of reverse search
                                   //since it may be that in front
                                   //of a gene there may be no tfbs
                                   //or no next gene
  while(rii !=jj)//begin not yet reached
    {
      if(IsTFBS(*rii))
	{
	  //counter++;
	  rii++;//continue the search
	}
      else
	jj=rii;//stops the search as rii==jj
               //makes sure that jj points to
               //something just past last found
               //tfbs (if any)
               //which may be either next gene
               //or begin of genome
    }
  return jj.base();//doing a base() operation on a reverse iterator
                   //gives a normal operator pointing to the element
                   //IN FRONT of the element pointed to by the reverse
                   //iterator;
}


Genome::iter Genome::FindRandomGenePosition() const
{
  std::list< iter > pos;
  std::list< iter >::iterator ipos;
  iter i;
  iter ii;  
  int randpos;
 
  if(gnrgenes_==0)
    {
      return (*ChromBBList).end();
    }
  else
    {
      i=(*ChromBBList).begin();
      while(i!=(*ChromBBList).end())
	{
	  if(IsGene(*i))
	    {
	      pos.push_back(i);
	    }
	  i++;
	}
      randpos=(int)(uniform()*gnrgenes_);
      if(randpos>gnrgenes_)
	{
	  printf("nonpresent gene selected\n");
	  exit(1);
	}
      ii=(*boost::next(pos.begin(),randpos));
      
      //ipos=pos.erase(pos.begin(),pos.end());
      return ii;
    }

}


Genome::iter Genome::GeneMutate(iter ii)
{
  double uu;
  iter first;
  iter last;
  iter jj;
  iter kk;
  int copylength;
  //double temp;
  Gene *gene;

  uu=uniform();
  if(uu<probnontandgenedupl)
    {
      last=ii;//gene position
      last++;//one further than the gene position
      first=FindFirstTFBSInFrontOfGene(ii);//first tfbs in front of gene
      copylength=distance(first,last);
      std::list<ChromBB*> ChromBBListTemp;//temporary list
      CopyPartOfGenomeToTempList(first,last,ChromBBListTemp);
      kk=FindRandomGenePosition();//find position of gene to insert in front of
      kk=FindFirstTFBSInFrontOfGene(kk);//find fist tfbs in front of this gene
      (*ChromBBList).splice(kk,ChromBBListTemp);//splice temporary list into chromosome
      glength_+=copylength;
      gnrtfbs_+=(copylength-1);
      gnrgenes_++;
      //make sure ii points to one position further than just duplicated gene
      ii=last;

    }
  else if(uu<probnontandgenedupl+probgenedel)
    {
      if(OnlyCopy(ii)==1)//only allowed to remove gene if there are multiple copies
	{
	  //delete gene and its upfront TFBS 
	  last=ii;//gene position
	  last++;//one further than the gene position
	  first=FindFirstTFBSInFrontOfGene(ii);//first tfbs in front of gene
	  jj=first;
	  while( jj != last ) 
	    {
	      gnrtfbs_--;
	      glength_--;
	      delete *jj;
	      jj++;
	    }
	  gnrgenes_--;//you know one gene is removed
	  gnrtfbs_++;//of the elements counted/removed above one was not tfbs but gene
	  ii=(*ChromBBList).erase(first,last);
	}
      else
	ii++;
    }
  else
    {
      ii++;
    }

  return ii;
}


int Genome::OnlyCopy(iter ii)
{
  Gene *gene;
  int type;
  Gene *gene2;
  int type2;
  iter jj;
  int secondcopy=0;
  
  gene=dynamic_cast<Gene *>(*ii);
  type=gene->type;
  
  for(jj=ChromBBList->begin();jj!=ChromBBList->end();jj++)
    {
      if(jj!=ii)//other position
	{
	  if(IsGene(*jj))//is it a gene
	    {
	      gene2=dynamic_cast<Gene *>(*jj);
	      type2=gene2->type;
	      
	      if(type2==type)
		{
		  secondcopy=1;
		  break;
		}
	    }
	}
    }

  return secondcopy;
}


Genome::iter Genome::TFBSMutate(iter ii)//OK
{
  double uu;
  TFBS *tfbs;
  TFBS *tfbsnew;
  //Gene *gene;
  //int genenr;
  iter jj;
  //double temp;
  //int inttemp;
  //int type;

  tfbs=dynamic_cast<TFBS *>(*ii);
  uu=uniform();
  if(uu<probtfbsweightrev)
    {
      (*tfbs).weight=-(*tfbs).weight;
      ii++;
    }
  else if(uu<probtfbsweightrev+probnontandtfbsdupl)
    {
      tfbsnew=new TFBS(*tfbs);
      jj=FindRandomGenePosition();
      jj=(*ChromBBList).insert(jj,tfbsnew);
      ii++;
      gnrtfbs_++;
      glength_++;
    }
  else if(uu<probtfbsweightrev+probnontandtfbsdupl+probtfbsdel)
    {
      delete (tfbs);
      ii=(*ChromBBList).erase(ii);
      gnrtfbs_--;
      glength_--;
    }
  else if(uu<probtfbsweightrev+probnontandtfbsdupl+probtfbsdel+probtfbstypeswitch)
    {
      (*tfbs).type=(int)(uniform()*NrGeneTypes);
      ii++;
    }
  else 
    ii++;

  return ii;
}


int Genome::CheckViability()
{
  int ok;
  iter i;
  int checkpresence[NrGeneTypes];
  int I;
  Gene *gene;
  int type;

  for(I=0;I<NrGeneTypes;I++)
    checkpresence[I]=0;

  ok=1;
  for(i=ChromBBList->begin();i!=ChromBBList->end();i++)
    {
      if(IsGene(*i))
	{
	  gene=dynamic_cast<Gene *>(*i);
	  type=gene->type;
	  checkpresence[type]++;
	}
    }
  for(I=0;I<NrGeneTypes;I++)
    {
      if(checkpresence[I]==0)
	{
	  ok=0;
	  break;
	}
    }
  if(ok==1)
    return TRUE;
  else
    return FALSE;
}
  


void Genome::ListContent()
{
  iter i;
  Gene *gene;
  TFBS *tfbs;

  printf("size list is %i\n",(*ChromBBList).size());

  for(i=(*ChromBBList).begin();i!=(*ChromBBList).end();i++)
     {
       if(IsGene(*i))
	 {
	   gene=dynamic_cast<Gene *>(*i);
	   printf("g%i",(*gene).type);
	 }
       else if(IsTFBS(*i))
	 {
	   tfbs=dynamic_cast<TFBS *>(*i);
	   if(tfbs->weight==-1)
	     printf("-%i",(*tfbs).type);
	   else
	     printf("+%i",(*tfbs).type);
	 }
     }
  printf("\n");
}

void Genome::PruneGenome()
{
  int nr;
  double uni;

  if(uniform()<0.4)//delete gene
    {
      uni=uniform();
      nr=(int)(uni*gnrgenes_);
      DeleteGene(nr);
      //printf("deleting gene nr %d, nr genes left:%d, nr tfbs left:%d \n", nr, gnrgenes_, gnrtfbs_);
      
    }
  else//delete tfbs
    { 
      uni=uniform();
      nr=(int)(uni*gnrtfbs_);
      DeleteTFBS(nr);
      //printf("deleting tfbs nr %d,  nr genes left:%d, nr tfbs left:%d\n", nr,gnrgenes_,gnrtfbs_);
    }
}

void Genome::RemoveSegmentationGene()
{
  iter it;
  iter last;
  iter first;
  iter jj;
  int counter=0;

  it=ChromBBList->begin();
  while(it!=ChromBBList->end())
    {
      if(IsGene(*it))
	{
	  if((*it)->type==SegmGeneNr)
	    {
	      //delete gene and its upfront TFBS 
	      last=it;//gene position
	      last++;//one further than the gene position
	      first=FindFirstTFBSInFrontOfGene(it);//first tfbs in front of gene
	      jj=first;
	      while( jj != last ) 
		{
		  gnrtfbs_--;
		  glength_--;
		  delete *jj;
		  jj++;
		}
	      gnrgenes_--;//you know one gene is removed
	      gnrtfbs_++;//of the elements counted/removed above one was not tfbs but gene
	      it=(*ChromBBList).erase(first,last);
	    }
	  else
	    counter++;
	}
      it++;
    }

}

void Genome::DeleteGene(int nr)
{
  iter it;
  iter last;
  iter first;
  iter jj;
  int counter=0;

  it=ChromBBList->begin();
  while(it!=ChromBBList->end())
    {
      if(IsGene(*it))
	{
	  if(counter==nr)
	    {
	      //delete gene and its upfront TFBS 
	      last=it;//gene position
	      last++;//one further than the gene position
	      first=FindFirstTFBSInFrontOfGene(it);//first tfbs in front of gene
	      jj=first;
	      while( jj != last ) 
		{
		  gnrtfbs_--;
		  glength_--;
		  delete *jj;
		  jj++;
		}
	      gnrgenes_--;//you know one gene is removed
	      gnrtfbs_++;//of the elements counted/removed above one was not tfbs but gene
	      it=(*ChromBBList).erase(first,last);

	      break;
	    }
	  else
	    counter++;
	}
      it++;
    }

}

void Genome::DeleteTFBS(int nr)
{
  iter it;
  int counter=0;

  it=ChromBBList->begin();
  while(it!=ChromBBList->end())
    {
      if(IsTFBS(*it))
	{
	  if(counter==nr)
	    {
	      delete (*it);
	      it=(*ChromBBList).erase(it);
	      gnrtfbs_--;
	      glength_--;

	      break;
	    }
	  else
	    counter++;
	}
      it++;
    }
}