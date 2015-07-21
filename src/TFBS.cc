#include "TFBS.hh"
#include "Header.hh"

TFBS::TFBS() : ChromBB() 
{
  type=0;
  weight=0;
  HH=0;
}

ChromBB* TFBS::clone() const
{
  return new TFBS(*this);
}

TFBS::~TFBS()
{
}

TFBS::TFBS(int t,int w,double H) : ChromBB() 
{
  type=t;
  weight=w;
   HH=H;
}

TFBS::TFBS(const TFBS &tfbs) : ChromBB(tfbs)
{
  type=tfbs.type;
  weight=tfbs.weight;
  HH=tfbs.HH;
  
}          
     
bool TFBS::isEqual( const ChromBB &a ) const
{
  const TFBS *tfbs =  dynamic_cast< const TFBS* >(& a );
  if( tfbs != NULL ) 
    return (this->type == tfbs->type && this->weight == tfbs->weight);
  else 
    return false;
}
