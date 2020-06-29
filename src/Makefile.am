bin_PROGRAMS=tnxverify
tnxverify_SOURCES=AMath.cpp ATimeSpace.cpp ACatalog.cpp ACatUCAC4.cpp ProjectTNX.cpp tnxverify.cpp

if DEBUG
  AM_CFLAGS = -g3 -O0 -Wall -DNDEBUG
  AM_CXXFLAGS = -g3 -O0 -Wall -DNDEBUG
else
  AM_CFLAGS = -O3 -Wall
  AM_CXXFLAGS = -O3 -Wall
endif

tnxverify_LDFLAGS = -L/usr/local/lib
BOOST_LIBS = -lboost_system-mt -lboost_filesystem-mt -lboost_date_time
tnxverify_LDADD = ${BOOST_LIBS} -lm -lcfitsio
