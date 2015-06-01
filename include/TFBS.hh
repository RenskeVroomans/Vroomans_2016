#ifndef TFBSHeader
#define TFBSHeader

#include "ChromBB.hh"

class TFBS : public ChromBB
{
 public:
  int weight;//+1 or -1
 
  TFBS();//constructor
  virtual ChromBB* clone() const;
  virtual ~TFBS();//destructor
  TFBS(int t,int w);//constructor with args
  explicit TFBS(const TFBS &tfbs);//copy constructor
  
  virtual bool isEqual( const ChromBB & ) const;
};
#endif
