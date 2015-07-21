#include "Gene.hh"
#include "Header.hh"

Gene::Gene() : ChromBB()
{
  type=0;
  DD=0;
  EE=0;
 
}

ChromBB* Gene::clone() const
{
  return new Gene(*this);
}

Gene::~Gene()
{
}

Gene::Gene(int t,double D, double E) : ChromBB()
{
  type=t;
  DD=D;
  EE=E;
}

Gene::Gene(const Gene &gene) : ChromBB(gene)
{
  type=gene.type;
  DD=gene.DD;
  EE=gene.EE;
}   
                     
bool Gene::isEqual( const ChromBB &a ) const
{
  const Gene *gene =  dynamic_cast< const Gene* >(& a );
  if( gene != NULL )
    return (this->type == gene->type );
  else 
    return false;
}