#ifndef ObjectStore_hh
#define ObjectStore_hh

#include <TObject.h>

#include <map>

class ObjectStore : public TObject
{
public:

  ObjectStore();

  virtual ~ObjectStore();

  static void clear();

protected:

private:

  void register_();
  void unregister_();

  static bool _clear;
  static std::map< TObject*, TObject* > _store;

};

#endif

