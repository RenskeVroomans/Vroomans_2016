#ifndef ChromBBHeader
#define ChromBBHeader
class ChromBB
{
 public:
  int type;

  ChromBB();
  virtual ChromBB *clone() const=0;
  virtual ~ChromBB();
  explicit ChromBB(const ChromBB &cbb);
  virtual bool isEqual( const ChromBB & ) const= 0;
};

inline bool operator==( const ChromBB &a, const ChromBB &b ) 
{ return a.isEqual( b ); }

inline bool operator!=( const ChromBB &a, const ChromBB &b ) 
{ return !a.isEqual( b ); }
#endif
