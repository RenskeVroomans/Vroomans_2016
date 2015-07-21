#ifndef GeneHeader
#define GeneHeader

#include "ChromBB.hh"

class Gene : public ChromBB
{
 public:
  double DD;
  double EE;
   
  Gene();//constructor
  virtual ChromBB* clone() const;
  ~Gene();//destructor
  Gene(int t,double D, double E);//constructor with arguments
  explicit Gene(const Gene &gene);//copy constructor
  
  virtual bool isEqual( const ChromBB & ) const;

};
#endif
