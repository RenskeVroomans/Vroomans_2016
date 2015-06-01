#ifndef GeneHeader
#define GeneHeader

#include "ChromBB.hh"

class Gene : public ChromBB
{
 public:
  
  Gene();//constructor
  virtual ChromBB* clone() const;
  ~Gene();//destructor
  Gene(int t);//constructor with arguments
  explicit Gene(const Gene &gene);//copy constructor
  
  virtual bool isEqual( const ChromBB & ) const;

};
#endif
