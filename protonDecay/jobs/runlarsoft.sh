#!/bin/sh
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v05_13_00 -q e9:prof
rundir=`pwd`
mkdir larDev
cd larDev
mrb newDev
cd srcs
source $rundir/larDev/localProducts_larsoft_v05_13_00_e9_prof/setup
mrbslp
mrb newProduct -c mulengthfilter
cd mulengthfilter/mulengthfilter
curl http://www.hep.shef.ac.uk/tmp/MuLengthFilter.tar.gz | tar xzf -
ls -lrt *
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j1
cd $rundir
cat > job.fcl << EOF
#include "services_dune.fcl"
#include "MUSUN.fcl"
#include "largeantmodules_dune.fcl"
#include "photpropservices_dune.fcl"
#include "detsimmodules_dune.fcl"
#include "opticaldetectormodules_dune.fcl"
#include "caldata_dune.fcl"
#include "hitfindermodules_dune.fcl"
#include "cluster_dune.fcl"
#include "trackfindermodules_dune.fcl"
#include "pandoramodules_dune.fcl"
#include "calorimetry_dune35t.fcl"
#include "mctrutht0matching.fcl"
#include "t0reco.fcl"
#include "opticaldetectormodules_dune.fcl"

process_name: MUSUNSim

services:
{
  # Load the service that manages root files for histograms.
  TFileService: { fileName: "MUSUNSim_hist.root" }
  TimeTracker:       {}
  Timing:       {}
  SimpleMemoryCheck:     { ignoreTotal: 1 } # default is one
  RandomNumberGenerator: {} #ART native random number generator
  @table::dunefd_simulation_services
  Geometry:     @local::dune10kt_geo
  message:              @local::dune_message_services_prod_debug
  FileCatalogMetadata:  @local::art_file_catalog_mc
                        @table::dunefd_services
}


#Start each new event with an empty event.
source:
{
  module_type: EmptyEvent
  timestampPlugin: { plugin_type: "GeneratedEventTimestamp" }
  maxEvents:   NEVENTS          # Number of events to create
  firstRun:    RUN           # Run number to use for this file
  firstEvent:  FIRST           # number of first event in the file
}

# Define and configure some modules to do work on each event.
# First modules are defined; they are scheduled later.
# Modules are grouped by type.
physics:
{

 producers:
 {
   rns:       { module_type: "RandomNumberSaver" }
   generator: @local::standard_MUSUN
   largeant:  @local::dunefd_largeant
   daq:       @local::dunefd_simwire
#   simcounter:     @local::dunefd_simcounter
   opdigi:         @local::dunefd_opdigi

 }

 filters:
 {
    mulen:      { module_type: "MuLengthFilter" }
 }
 #define the producer and filter modules for this path, order matters, 
 #filters reject all following items.  see lines starting physics.producers below
 simulate: [ rns, generator, largeant, mulen
             ]
 
 #define the output stream, there could be more than one if using filters 
 stream1:  [ out1 ]
 #stream1:  [ ]

 #trigger_paths is a keyword and contains the paths that modify the art::event, 
 #ie filters and producers
 trigger_paths: [simulate] 

 #end_paths is a keyword and contains the paths that do not modify the art::Event, 
 #ie analyzers and output streams.  these all run simultaneously
 end_paths:     [stream1]  
}


#block to define where the output goes.  if you defined a filter in the physics
#block and put it in the trigger_paths then you need to put a SelectEvents: {SelectEvents: [XXX]}
#entry in the output stream you want those to go to, where XXX is the label of the filter module(s)
outputs:
{
 out1:
 {
   module_type: RootOutput
   fileName:    "nodaq.root"
   dataTier:    "simulation"
   compressionLevel: 1
   SelectEvents: { SelectEvents: [ simulate ] }
 }
}
EOF
lar -c job.fcl
