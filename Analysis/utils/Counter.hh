#ifndef Counter_hh
#define Counter_hh

//
// original authors: Bob Jacobsen and Gautier Hamel de Monchenault
//   (Lawrence Berkeley Laboratory)

template <class T>
class Counter 
{
public:

  static int nAlloc() { return _nAlloc; }

  Counter()    { _nAlloc++; } 
  ~Counter()   { --_nAlloc; }

private:

  static int _nAlloc;

  int dum() const;
};

#endif

