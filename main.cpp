#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <iostream>
#include <mutex>
#include <osg/ArgumentParser>
#include <osg/Array>
#include <osg/Callback>
#include <osg/PositionAttitudeTransform>
#include <osg/Texture2D>
#include <osgGA/GUIEventHandler>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <thread>
#include <vector>

#define dimof(x) (sizeof (x) / sizeof (*x))

enum sort_id_t
{
  SID_UNSPECIFIED,
  SID_INSERTION,
  SID_BUBBLE,
  SID_COMB,
  SID_SHELL,
  SID_SELECTION,
  SID_QUICK
};

namespace
{

struct meta_t
{
  int key_;
  sort_id_t id_;
  const char *name_;
} static const g_meta[] = {
  { osgGA::GUIEventAdapter::KEY_1, SID_INSERTION, "Insertion sort" },
  { osgGA::GUIEventAdapter::KEY_2, SID_BUBBLE, "Bubble sort" },
  { osgGA::GUIEventAdapter::KEY_3, SID_COMB, "Comb sort" },
  { osgGA::GUIEventAdapter::KEY_4, SID_SHELL, "Shell sort" },
  { osgGA::GUIEventAdapter::KEY_5, SID_SELECTION, "Selection sort" },
  { osgGA::GUIEventAdapter::KEY_6, SID_QUICK, "Quicksort" },
};

} // anonymous namespace

static void
modify_image (osg::Image *img, std::vector <int> *a)
{
  // ABGR.
  static const uint32_t colors[] = {
    0xFFFF0000,
    0xFF800000,
    0xFF400000,
    0xFFFFFF00,
    0xFF808000,
    0xFF404000,
    0xFF00FF00,
    0xFF008000,
    0xFF004000,
    0xFF00FFFF,
    0xFF008080,
    0xFF004040,
    0xFFFF00FF,
    0xFF800080,
    0xFF400040,
    0xFF0000FF,
    0xFF000080,
    0xFF000040,
  };

  int width = static_cast <int> (a->size ());
  uint32_t *data = reinterpret_cast <uint32_t *> (img->data ());
  for (auto &e : *a)
  {
    double f = (e - 1.0) / (width - 1.0) * dimof (colors);
    size_t ci = std::min (dimof (colors) - 1, static_cast <size_t> (f));
    uint32_t color = colors[ci];
    std::fill (data, data + e, color);
    std::fill (data + e, data + width, 0x00FFFFFF);
    data += width;
  }
}

static osg::Image *
create_image (std::vector <int> *a)
{
  int image_size = static_cast <int> (a->size ());
  osg::ref_ptr <osg::Image> img (new osg::Image ());
  img->allocateImage (image_size, image_size, 1, GL_RGBA, GL_UNSIGNED_BYTE);
  return img.release ();
}

namespace
{

class context_t
{
public:
  explicit
  context_t (int sort_count);

  std::mutex m_;
  std::condition_variable cv_;
  std::condition_variable reset_cv_;
  bool stopped_;
  bool reset_;
  size_t update_counter_;
  size_t render_counter_;
  std::vector <int> a_;
  sort_id_t id_;

private:
  context_t (const context_t &other) = delete;

  context_t &
  operator = (const context_t &other) = delete;
};

} // anonymous namespace

context_t::context_t (int sort_count)
: m_ (),
  cv_ (),
  reset_cv_ (),
  stopped_ (false),
  reset_ (false),
  update_counter_ (0),
  render_counter_ (1),
  a_ (),
  id_ (SID_UNSPECIFIED)
{
  if (sort_count < 0)
    throw std::invalid_argument (__FUNCTION__);
  a_.resize (static_cast <size_t> (sort_count));
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
    for (size_t i = 0; i < dimof (g_meta); ++i)
    {
      if (ea.getKey () == g_meta[i].key_)
      {
        std::unique_lock <std::mutex> lock (ctx_->m_);
        ctx_->reset_ = true;
        ctx_->id_ = g_meta[i].id_;
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
  modify_image (img, &ctx_->a_);
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
create_geom (float aspect_ratio, float half_width)
{
  // Oriented as though we are looking towards the positive Y-axis where
  // Positive Z is up, and positive X is right.

  float half_height = half_width * aspect_ratio;

  osg::ref_ptr <osg::Vec3Array> v (new osg::Vec3Array ());
  v->push_back (osg::Vec3 (-half_width, 0.0f, -half_height));
  v->push_back (osg::Vec3 (half_width, 0.0f, -half_height));
  v->push_back (osg::Vec3 (half_width, 0.0f, half_height));
  v->push_back (osg::Vec3 (-half_width, 0.0f, half_height));

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
create_model (context_t *ctx, float aspect_ratio, float half_width)
{
  osg::ref_ptr <osg::Image> img (create_image (&ctx->a_));
  osg::ref_ptr <osg::Texture2D> texture (create_texture ());
  osg::ref_ptr <texture_cb_t> tcb (new texture_cb_t (ctx));
  texture->setUpdateCallback (tcb);
  texture->setDataVariance (osg::Object::DYNAMIC);
  texture->setImage (img);

  osg::ref_ptr <osg::Geometry> geom (create_geom (aspect_ratio, half_width));
  osg::StateSet *ss = geom->getOrCreateStateSet ();
  ss->setTextureAttributeAndModes (0, texture, osg::StateAttribute::ON); 
  ss->setMode (GL_BLEND, osg::StateAttribute::ON);
  ss->setMode (GL_LIGHTING, osg::StateAttribute::OFF);

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

#define WAIT_FOR_MODEL_UPDATE(ctx, lock) do { \
  ctx->cv_.wait_for (lock, std::chrono::milliseconds (500), [ctx] { \
    return ctx->render_counter_ == ctx->update_counter_; \
  }); \
  ++ctx->update_counter_; \
} while (0);


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
      WAIT_FOR_MODEL_UPDATE (ctx, lock);
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
        WAIT_FOR_MODEL_UPDATE (ctx, lock);
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
          WAIT_FOR_MODEL_UPDATE (ctx, lock)
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
shellsort (IT lo, IT hi, context_t *ctx)
{
  static const int gaps[] = { 32003, 1701, 701, 301, 132, 57, 23, 10, 4, 1 };

  for (auto &gap : gaps)
  {
    if (gap >= std::distance (lo, hi))
      continue;

    for (IT i = lo + gap; i != hi; ++i)
    {
      auto current = *i;
      IT j = i;
      while (std::distance (lo, j-gap) >= 0 && current < *(j-gap))
      {
        std::unique_lock <std::mutex> lock (ctx->m_);
        if (ctx->reset_)
          return;
        WAIT_FOR_MODEL_UPDATE (ctx, lock)
        *j = *(j-gap);
        j -= gap;
      }
      std::unique_lock <std::mutex> lock (ctx->m_);
      if (ctx->reset_)
        return;
      WAIT_FOR_MODEL_UPDATE (ctx, lock)
      *j = current;
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
      WAIT_FOR_MODEL_UPDATE (ctx, lock)
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
    if (ctx->reset_)
      return lo;
    WAIT_FOR_MODEL_UPDATE (ctx, lock)
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
      for (size_t i = 0; i < ctx->a_.size (); ++i)
        ctx->a_[i] = (rand () % ctx->a_.size ()) + 1;
      ctx->reset_ = false;
      ctx->update_counter_ = 0;
      id = ctx->id_;
    }

    std::this_thread::sleep_for (std::chrono::milliseconds (100));
    switch (id)
    {
    case SID_INSERTION:
      insertionsort (begin (ctx->a_), end (ctx->a_), ctx);
      break;
    case SID_BUBBLE:
      bubblesort (begin (ctx->a_), end (ctx->a_), ctx);
      break;
    case SID_COMB:
      combsort (begin (ctx->a_), end (ctx->a_), ctx);
      break;
    case SID_SHELL:
      shellsort (begin (ctx->a_), end (ctx->a_), ctx);
      break;
    case SID_SELECTION:
      selectionsort (begin (ctx->a_), end (ctx->a_), ctx);
      break;
    case SID_QUICK:
      quicksort2 (begin (ctx->a_), end (ctx->a_), ctx);
      break;
    case SID_UNSPECIFIED:
    default:
      {
        std::unique_lock <std::mutex> lock (ctx->m_);
        std::sort (begin (ctx->a_), end (ctx->a_));
      }
      break;
    }

    std::unique_lock <std::mutex> lock (ctx->m_);
    ctx->reset_cv_.wait (lock, [ctx] { return ctx->reset_; });
  }
}

static float
get_aspect_ratio (osgViewer::Viewer &viewer)
{
  osgViewer::ViewerBase::Windows windows;
  viewer.getWindows (windows);
  int width = 1920;
  int height = 1080;
  if (!windows.empty ())
  {
    int x, y;
    windows[0]->getWindowRectangle (x, y, width, height);
  }
  return 1.0f * height / width;
}

static bool
is_valid_sort_count (int sort_count)
{
  for (size_t i = 1; i <= 12; ++i)
  {
    if (sort_count == pow (2, i))
      return true;
  }
  return false;
}

static void
setup_usage (osg::ApplicationUsage *usage, const std::string &appname)
{
  usage->setApplicationName (appname);
  usage->setDescription (appname + " visualizes sorting algorithms");
  usage->addCommandLineOption("--count", "How many numbers to sort", "128");
  for (size_t i = 0; i < dimof (g_meta); ++i)
    usage->addKeyboardMouseBinding (g_meta[i].key_, g_meta[i].name_);

  usage->addKeyboardMouseBinding (osgGA::GUIEventAdapter::KEY_R, "Restart Sort");
}

int
main (int argc, char *argv[])
{
  static const float HALF_WIDTH = 100.0f;
  static const double HALF_VIEW_WIDTH = HALF_WIDTH * 1.05;

  osg::ArgumentParser args (&argc, argv);
  setup_usage (args.getApplicationUsage (), args.getApplicationName ());

  osgViewer::Viewer viewer (args);

  unsigned int helpType;
  if ((helpType = args.readHelpType ()))
  {
    args.getApplicationUsage ()->write (std::cout, helpType);
    return 2;
  }

  int sort_count = 128;
  args.read("--count", sort_count);
  if (!is_valid_sort_count (sort_count))
    args.reportError ("Invalid count: try powers of 2 up to 4096");

  args.reportRemainingOptionsAsUnrecognized ();
  if (args.errors ())
  {
    args.writeErrorMessages (std::cout);
    return 2;
  }

  viewer.realize ();

  // Aspect ratio as a fraction.
  float aspect_ratio = get_aspect_ratio (viewer);

  context_t ctx (sort_count);

  osg::ref_ptr <osg::Group> model = create_model (
    &ctx, aspect_ratio, HALF_WIDTH
  );
  viewer.setSceneData (model);

  osg::ref_ptr <event_handler_t> eh (new event_handler_t (&ctx));
  viewer.addEventHandler (eh);
  viewer.addEventHandler (new osgViewer::StatsHandler ());

  auto cam = viewer.getCamera ();
  cam->setClearColor (osg::Vec4 (0.1f, 0.1f, 0.1f, 1.0f));
  cam->setProjectionMatrixAsOrtho2D (
    -HALF_VIEW_WIDTH, HALF_VIEW_WIDTH,
    -HALF_VIEW_WIDTH * aspect_ratio, HALF_VIEW_WIDTH * aspect_ratio
  );

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
