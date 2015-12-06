#include <cstdio>
#include <osg/ArgumentParser>
#include <osg/Array>
#include <osg/Callback>
#include <osg/PositionAttitudeTransform>
#include <osg/Texture2D>
#include <osgGA/GUIEventHandler>
#include <osgViewer/Viewer>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#define dimof(x) (sizeof (x) / sizeof (*x))

enum sort_id_t
{
  SID_UNSPECIFIED,
  SID_INSERTION,
  SID_BUBBLE,
  SID_COMB,
  SID_SELECTION,
  SID_QUICK
};

static void
modify_image (osg::Image *img, std::vector <int> *a)
{
  static const int IMAGE_WIDTH = 256;
  static const int IMAGE_HEIGHT = 256;

  // ABGR.
  static const uint32_t colors[] = {
    0xFFFF0000,
    0xFF800000,
    0xFFFFFF00,
    0xFF808000,
    0xFF00FF00,
    0xFF00FFFF,
    0xFF008080,
    0xFF008000,
    0xFF0000FF,
    0xFF000080,
  };

  uint32_t *data = reinterpret_cast <uint32_t *> (img->data ());
  for (int i = 0; i < IMAGE_HEIGHT; ++i)
  {
    for (int j = 0; j < IMAGE_WIDTH; ++j)
    {
      uint32_t color = 0xFF000000;
      if (i < a->size ())
      {
        if (j <= (*a)[i])
          color = colors[size_t (0.5 + (((*a)[i] / 256.0) * (dimof (colors) - 1)))];
      }
      size_t index = IMAGE_WIDTH * i + j;
      data[index] = color;
    }
  }
}

static osg::Image *
create_image (std::vector <int> *a)
{
  static const int IMAGE_WIDTH = 256;
  static const int IMAGE_HEIGHT = 256;
  osg::ref_ptr <osg::Image> img (new osg::Image ());
  img->allocateImage (IMAGE_WIDTH, IMAGE_HEIGHT, 1, GL_RGBA, GL_UNSIGNED_BYTE);
  return img.release ();
}

namespace
{

class context_t
{
public:
  context_t (std::vector <int> *a);

  std::mutex m_;
  std::condition_variable cv_;
  std::condition_variable reset_cv_;
  bool stopped_;
  bool reset_;
  size_t update_counter_;
  size_t render_counter_;
  std::vector <int> *a_;
  sort_id_t id_;

private:
  context_t (const context_t &other) = delete;

  context_t &
  operator = (const context_t &other) = delete;
};

} // anonymous namespace

context_t::context_t (std::vector <int> *a)
: m_ (),
  cv_ (),
  reset_cv_ (),
  stopped_ (false),
  reset_ (false),
  update_counter_ (0),
  render_counter_ (1),
  a_ (a),
  id_ (SID_UNSPECIFIED)
{
}

namespace
{

class event_handler_t : public osgGA::GUIEventHandler
{
public:
  explicit
  event_handler_t (context_t *ctx);

  virtual bool
  handle (
    const osgGA::GUIEventAdapter &ea, osgGA::GUIActionAdapter &aa,
    osg::Object *obj, osg::NodeVisitor *nv
  );

private:
  context_t *ctx_;
};

} // anonymous namespace

event_handler_t::event_handler_t (context_t *ctx)
: osgGA::GUIEventHandler (),
  ctx_ (ctx)
{
}

bool
event_handler_t::handle (
  const osgGA::GUIEventAdapter &ea, osgGA::GUIActionAdapter &aa,
  osg::Object *obj, osg::NodeVisitor *nv )
{
  struct map_t
  {
    int key_;
    sort_id_t id_;
  } static const m[] = {
    { osgGA::GUIEventAdapter::KEY_1, SID_INSERTION },
    { osgGA::GUIEventAdapter::KEY_2, SID_BUBBLE },
    { osgGA::GUIEventAdapter::KEY_3, SID_COMB },
    { osgGA::GUIEventAdapter::KEY_4, SID_SELECTION },
    { osgGA::GUIEventAdapter::KEY_5, SID_QUICK },
  };

  if (osgGA::GUIEventAdapter::KEYDOWN == ea.getEventType ())
  {
    if (osgGA::GUIEventAdapter::KEY_R == ea.getKey ())
    {
      std::unique_lock <std::mutex> lock (ctx_->m_);
      ctx_->reset_ = true;
      lock.unlock ();
      ctx_->reset_cv_.notify_one ();
      return true;
    }
    for (size_t i = 0; i < dimof (m); ++i)
    {
      if (ea.getKey () == m[i].key_)
      {
        std::unique_lock <std::mutex> lock (ctx_->m_);
        ctx_->reset_ = true;
        ctx_->id_ = m[i].id_;
        lock.unlock ();
        ctx_->reset_cv_.notify_one ();
        return true;
      }
    }
  }

  return false;
}

namespace
{

class texture_cb_t : public osg::StateAttributeCallback
{
public:
  explicit
  texture_cb_t (context_t *ctx);

  virtual
  ~texture_cb_t ();

  virtual void
  operator () (osg::StateAttribute *sa, osg::NodeVisitor *nv);

private:
  context_t *ctx_;
};

} // anonymous namespace

texture_cb_t::texture_cb_t (context_t *ctx)
: Callback (),
  ctx_ (ctx)
{
}

texture_cb_t::~texture_cb_t ()
{
}

void
texture_cb_t::operator () (osg::StateAttribute *sa, osg::NodeVisitor *nv)
{
  if (!sa)
    return;

  osg::Texture *tex = sa->asTexture ();
  if (!tex)
    return;

  if (GL_TEXTURE_2D != tex->getTextureTarget ())
    return;

  std::unique_lock <std::mutex> lock (ctx_->m_);
  if (ctx_->stopped_ || ctx_->render_counter_ == ctx_->update_counter_)
  {
    lock.unlock ();
    ctx_->cv_.notify_one ();
    return;
  }

  osg::Image *img = tex->getImage (0);
  modify_image (img, ctx_->a_);
  img->dirty ();

  ctx_->render_counter_ = ctx_->update_counter_;
  lock.unlock ();
  ctx_->cv_.notify_one ();
}

static osg::Geometry *
create_compass ()
{
  osg::ref_ptr <osg::Vec3Array> v (new osg::Vec3Array ());
  v->push_back (osg::Vec3 (0.0f, 0.0f, 0.0f));
  v->push_back (osg::Vec3 (10.0f, 0.0f, 0.0f));
  v->push_back (osg::Vec3 (0.0f, 0.0f, 0.0f));
  v->push_back (osg::Vec3 (0.0f, 10.0f, 0.0f));
  v->push_back (osg::Vec3 (0.0f, 0.0f, 0.0f));
  v->push_back (osg::Vec3 (0.0f, 0.0f, 10.0f));

  static const osg::Vec4 color (1.0f, 1.0f, 1.0f, 1.0f);
  osg::ref_ptr <osg::Vec4Array> c (new osg::Vec4Array ());
  c->reserve (v->size ());
  c->push_back (osg::Vec4 (0.0f, 0.0f, 1.0f, 1.0f));
  c->push_back (osg::Vec4 (0.0f, 0.0f, 1.0f, 1.0f));
  c->push_back (osg::Vec4 (0.0f, 1.0f, 0.0f, 1.0f));
  c->push_back (osg::Vec4 (0.0f, 1.0f, 0.0f, 1.0f));
  c->push_back (osg::Vec4 (1.0f, 0.0f, 0.0f, 1.0f));
  c->push_back (osg::Vec4 (1.0f, 0.0f, 0.0f, 1.0f));

  static const osg::Vec3 normal (0.0f, -1.0f, 0.0f);
  osg::ref_ptr <osg::Vec3Array> n (new osg::Vec3Array ());
  n->reserve (v->size ());
  for (size_t i = 0; i < v->size (); ++i)
    n->push_back (normal);

  osg::ref_ptr <osg::Geometry> out (new osg::Geometry ());
  out->setVertexArray (v);
  out->setColorArray (c, osg::Array::BIND_PER_VERTEX);
  out->setNormalArray (n, osg::Array::BIND_PER_VERTEX);
  out->addPrimitiveSet (new osg::DrawArrays (osg::PrimitiveSet::LINES, 0, 6));

  return out.release ();
}

static osg::Geometry *
create_geom ()
{
  // Oriented as though we are looking towards the positive Y-axis where
  // Positive Z is up, and positive X is right.
  osg::ref_ptr <osg::Vec3Array> v (new osg::Vec3Array ());
  v->push_back (osg::Vec3 (-100.0f, 0.0f, -64.0f));
  v->push_back (osg::Vec3 (100.0f, 0.0f, -64.0f));
  v->push_back (osg::Vec3 (100.0f, 0.0f, 64.0f));
  v->push_back (osg::Vec3 (-100.0f, 0.0f, 64.0f));

  static const osg::Vec4 color (1.0f, 1.0f, 1.0f, 1.0f);
  osg::ref_ptr <osg::Vec4Array> c (new osg::Vec4Array ());
  c->reserve (v->size ());
  for (size_t i = 0; i < v->size (); ++i)
    c->push_back (color);

  static const osg::Vec3 normal (0.0f, -1.0f, 0.0f);
  osg::ref_ptr <osg::Vec3Array> n (new osg::Vec3Array ());
  n->reserve (v->size ());
  for (size_t i = 0; i < v->size (); ++i)
    n->push_back (normal);

  osg::ref_ptr <osg::Vec2Array> tc (new osg::Vec2Array ());
  tc->reserve (v->size ());
  tc->push_back (osg::Vec2 (0.0f, 0.0f));
  tc->push_back (osg::Vec2 (0.0f, 1.0f));
  tc->push_back (osg::Vec2 (1.0f, 1.0f));
  tc->push_back (osg::Vec2 (1.0f, 0.0f));

  osg::ref_ptr <osg::Geometry> out (new osg::Geometry ());
  out->setVertexArray (v);
  out->setColorArray (c, osg::Array::BIND_PER_VERTEX);
  out->setNormalArray (n, osg::Array::BIND_PER_VERTEX);
  out->setTexCoordArray (0, tc, osg::Array::BIND_PER_VERTEX);
  out->addPrimitiveSet (new osg::DrawArrays (osg::PrimitiveSet::QUADS, 0, 4));

  return out.release ();
}

static osg::Texture2D *
create_texture ()
{
  osg::ref_ptr <osg::Texture2D> t (new osg::Texture2D ());
  t->setDataVariance (osg::Object::DYNAMIC);
  t->setFilter (osg::Texture::MIN_FILTER, osg::Texture::NEAREST);
  t->setFilter (osg::Texture::MAG_FILTER, osg::Texture::NEAREST);
  t->setWrap (osg::Texture::WRAP_S, osg::Texture::CLAMP);
  t->setWrap (osg::Texture::WRAP_T, osg::Texture::CLAMP);

  return t.release ();
}

static osg::Group *
create_model (context_t *ctx)
{
  osg::ref_ptr <osg::Image> img (create_image (ctx->a_));
  osg::ref_ptr <osg::Texture2D> texture (create_texture ());
  osg::ref_ptr <texture_cb_t> tcb (new texture_cb_t (ctx));
  texture->setUpdateCallback (tcb);
  texture->setDataVariance (osg::Object::DYNAMIC);
  texture->setImage (img);

  osg::ref_ptr <osg::Geometry> geom (create_geom ());
  osg::StateSet *ss = geom->getOrCreateStateSet ();
  ss->setTextureAttributeAndModes (0, texture, osg::StateAttribute::ON); 

  osg::ref_ptr <osg::Group> out (new osg::Group ());
  out->addChild (geom);

#if 0
  osg::ref_ptr <osg::PositionAttitudeTransform> pat (new osg::PositionAttitudeTransform ());
  pat->addChild (create_compass ());
  pat->setPosition (osg::Vec3 (-120.0f, 0.0f, 0.0f));
  out->addChild (pat);
#endif

  return out.release ();
}

template <typename IT> void
insertionsort (IT lo, IT hi, context_t *ctx)
{
  if (lo == hi)
    return;

  for (IT i = lo+1; i != hi; ++i)
  {
    for (IT j = i; j != lo && *j < *(j-1); --j)
    {
      std::unique_lock <std::mutex> lock (ctx->m_);
      if (ctx->reset_)
        return;
      ctx->cv_.wait_for (lock, std::chrono::milliseconds (10), [ctx] {
        return ctx->render_counter_ == ctx->update_counter_;
      });
      ++ctx->update_counter_;
      std::iter_swap (j, j-1);
    }
  }
}

template <typename IT> void
bubblesort (IT lo, IT hi, context_t *ctx)
{
  for (IT j = hi; lo != j; --j)
  {
    bool swapped = false;
    for (IT i = lo + 1; i != j; ++i)
    {
      if (*i < *(i-1))
      {
        std::unique_lock <std::mutex> lock (ctx->m_);
        if (ctx->reset_)
          return;
        ctx->cv_.wait_for (lock, std::chrono::milliseconds (10), [ctx] {
          return ctx->render_counter_ == ctx->update_counter_;
        });
        ++ctx->update_counter_;
        std::iter_swap (i, i-1);
        swapped = true;
      }
    }
    if (!swapped)
      break;
  }
}

template <typename IT> void
combsort (IT lo, IT hi, context_t *ctx)
{
  bool swapped = false;
  auto gap = std::distance (lo, hi);

  while (gap != 1 && !swapped)
  {
    gap = std::max (int (gap / 1.3), 1);
    for (IT j = hi; lo != j; --j)
    {
      bool swapped = false;
      for (IT i = lo; std::distance (i+gap,j) > 0; ++i)
      {
        if (*(i+gap) < *i)
        {
          std::unique_lock <std::mutex> lock (ctx->m_);
          if (ctx->reset_)
            return;
          ctx->cv_.wait_for (lock, std::chrono::milliseconds (100), [ctx] {
            return ctx->render_counter_ == ctx->update_counter_;
          });
          ++ctx->update_counter_;
          std::iter_swap (i, i+gap);
          swapped = true;
        }
      }
      if (!swapped)
        break;
    }
  }
}

template <typename IT> void
selectionsort (IT lo, IT hi, context_t *ctx)
{
  if (lo == hi)
    return;

  for (IT j = lo; j != hi-1; ++j)
  {
    std::this_thread::sleep_for (std::chrono::milliseconds (100));
    IT min = j;
    for (IT i = j+1; i != hi; ++i)
    {
      if (*i < *min)
        min = i;
    }
    if (j != min)
    {
      std::unique_lock <std::mutex> lock (ctx->m_);
      if (ctx->reset_)
        return;
      ctx->cv_.wait_for (lock, std::chrono::milliseconds (10), [ctx] {
        return ctx->render_counter_ == ctx->update_counter_;
      });
      ++ctx->update_counter_;
      std::iter_swap (j, min);
    }
  }
}

template <typename IT, typename P> IT
partition2 (IT lo, IT hi, P p, context_t *ctx)
{
  while (lo != hi)
  {
    while (p (*lo))
    {
      if (++lo == hi)
        return lo;
    }
    do
    {
      if (lo == --hi)
        return lo;
    } while (!p (*hi));
    std::unique_lock <std::mutex> lock (ctx->m_);
    ctx->cv_.wait_for (lock, std::chrono::milliseconds (100), [ctx] {
      return ctx->render_counter_ == ctx->update_counter_;
    });
    ++ctx->update_counter_;
    std::iter_swap (lo++, hi);
  }
  return lo;
}

template <typename IT> void
quicksort2 (IT lo, IT hi, context_t *ctx)
{
  if (lo == hi)
    return;

  {
    std::unique_lock <std::mutex> lock (ctx->m_);
    if (ctx->reset_)
      return;
  }

  auto t = *(lo + std::distance (lo, hi)/2);
  IT m1 = partition2 (lo, hi, [&t] (const auto &e) { return e < t; }, ctx);
  IT m2 = partition2 (m1, hi, [&t] (const auto &e) { return !(t < e); }, ctx);
  quicksort2 (lo, m1, ctx);
  quicksort2 (m2, hi, ctx);
}

static int
get_rand ()
{
  return (rand () % 256) + 1;
}

static void
sort (context_t *ctx)
{
  for (;;)
  {
    sort_id_t id = SID_UNSPECIFIED;

    {
      std::unique_lock <std::mutex> lock (ctx->m_);
      if (ctx->stopped_)
        return;
      for (size_t i = 0; i < ctx->a_->size (); ++i)
        (*ctx->a_)[i] = get_rand ();
      ctx->reset_ = false;
      ctx->update_counter_ = 0;
      id = ctx->id_;
    }

    std::this_thread::sleep_for (std::chrono::seconds (1));
    switch (id)
    {
    case SID_INSERTION:
      insertionsort (begin (*(ctx->a_)), end (*(ctx->a_)), ctx);
      break;
    case SID_BUBBLE:
      bubblesort (begin (*(ctx->a_)), end (*(ctx->a_)), ctx);
      break;
    case SID_COMB:
      combsort (begin (*(ctx->a_)), end (*(ctx->a_)), ctx);
      break;
    case SID_SELECTION:
      selectionsort (begin (*(ctx->a_)), end (*(ctx->a_)), ctx);
      break;
    case SID_QUICK:
      quicksort2 (begin (*(ctx->a_)), end (*(ctx->a_)), ctx);
      break;
    case SID_UNSPECIFIED:
    default:
      break;
    }

    std::unique_lock <std::mutex> lock (ctx->m_);
    while (!ctx->reset_)
    {
      ctx->reset_cv_.wait_for (
        lock, std::chrono::milliseconds (100), [ctx] { return ctx->reset_; }
      );
    }
  }
}

int
main (int argc, char *argv[])
{
  osg::ArgumentParser arguments (&argc, argv);
  osgViewer::Viewer viewer (arguments);

  std::vector <int> a (256);
  context_t ctx (&a);

  osg::ref_ptr <osg::Group> model = create_model (&ctx);
  viewer.setSceneData (model);

  osg::ref_ptr <event_handler_t> eh (new event_handler_t (&ctx));
  viewer.addEventHandler (eh);

  std::thread th (sort, &ctx);

  viewer.run ();

  {
    std::unique_lock <std::mutex> lock (ctx.m_);
    ctx.stopped_ = true;
    ctx.reset_ = true;
  }
  ctx.reset_cv_.notify_one ();

  th.join ();
  return 0;
}
