#include "Analysis/utils/ObjectStore.hh"

#include <cassert>

using namespace std;

bool ObjectStore::_clear = false;
map< TObject*, TObject* > ObjectStore::_store;

ObjectStore::ObjectStore()
{
  register_();
}

ObjectStore::~ObjectStore()
{
  unregister_();
}

void
ObjectStore::register_()
{
  assert( _store.count(this)==0 );
  _store[this] = this;
}

void
ObjectStore::unregister_()
{
  if( !_clear )
    {
      assert( _store.count(this)==1 );
      _store[this] = 0;
    }
}

void
ObjectStore::clear()
{
  _clear = true;
  for( map<TObject*,TObject*>::iterator it = _store.begin();
       it!= _store.end(); ++it )
    {
      delete it->second;
    }
  _store.clear();
  _clear = false;
}
