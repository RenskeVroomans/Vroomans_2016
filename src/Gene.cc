#include "Gene.hh"
#include "Header.hh"

Gene::Gene() : ChromBB()
{
  type=0;
}

ChromBB* Gene::clone() const
{
  return new Gene(*this);
}

Gene::~Gene()
{
}

Gene::Gene(int t) : ChromBB()
{
  type=t;
}

Gene::Gene(const Gene &gene) : ChromBB(gene)
{
  type=gene.type;
}   
                     
bool Gene::isEqual( const ChromBB &a ) const
{
  const Gene *gene =  dynamic_cast< const Gene* >(& a );
  if( gene != NULL )
    return (this->type == gene->type );
  else 
    return false;
}