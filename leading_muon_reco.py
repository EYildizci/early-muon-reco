from icecube import icetray, dataclasses, dataio, simclasses, MuonGun
import numpy as np
import matplotlib.path as mpltPath
from argparse import ArgumentParser
from icecube.icetray import I3Units
from I3Tray import I3Tray
import pandas as pd
import math
import copy


gcdFile = dataio.I3File('/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz')


frame = gcdFile.pop_frame()
while not frame.Has('I3Geometry'):
    frame = gcdFile.pop_frame()
geometry = frame.Get('I3Geometry')

om_geometry = geometry.omgeo


dom_x_positions=np.zeros((87,67))
dom_y_positions=np.zeros((87,67))
dom_z_positions=np.zeros((87,67))


for om, geo_info in om_geometry:
    dom_x_positions[om[0],om[1]]=geo_info.position.x
    dom_y_positions[om[0],om[1]]=geo_info.position.y
    dom_z_positions[om[0],om[1]]=geo_info.position.z
    
energy_range=np.loadtxt('./muon_range_1_to_3_w_150_steps')

def find_nearest(array,value):
    energies = array[:,0]
    ranges = array[:,1]
    idx = np.searchsorted(ranges, value, side="left")
    if idx > 0 and (idx == len(ranges) or math.fabs(value - ranges[idx-1]) < math.fabs(value - ranges[idx])):
        return energies[idx-1]
    else:
        return energies[idx]


dom_z_positions[0,:]=-2000



parser = ArgumentParser()
parser.add_argument("-i", "--infile", dest="infile",
                    default=None, type=str, required=True, 
                    help="the full path to the input file")

args = parser.parse_args()
infileNumber = args.infile

infile = '/data/sim/IceCube/2020/filtered/level2/neutrino-generator/21871/p0=0.0_p1=0.0_domeff=1.00/0000000-0000999/Level2_IC86.2020_NuE.021871.'+infileNumber+'.i3.zst'

tray = I3Tray()
tray.Add('I3Reader', Filename=infile)


def EarlyHitIdentifier(frame):
    event_id = frame['I3EventHeader'].event_id
    print(event_id)
    mcTree = frame['I3MCTree']
    #reco_early_muon_mcTree=copy.copy(mcTree)
    #reco_early_muon_mcTree.clear()
    mc_reco_early_muon_mcTree=copy.copy(mcTree)
    mc_reco_early_muon_mcTree.clear()
    #no_dir_mc_reco_early_muon_mcTree=copy.copy(mcTree)
    #no_dir_mc_reco_early_muon_mcTree.clear()

    #Create a list of muon and child ids
    muonIds = np.empty((0,1), int)
    muonChildIds = np.empty((0,1), int)

    try:
        hits = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'SplitInIcePulses')
    except:
        return False
    #reco_early_pulse_mask = dataclasses.I3RecoPulseSeriesMapMask(frame, 'SplitInIcePulses')
    mc_early_pulse_mask = dataclasses.I3RecoPulseSeriesMapMask(frame, 'SplitInIcePulses')
    mc_muons_pulse_mask = dataclasses.I3RecoPulseSeriesMapMask(frame, 'SplitInIcePulses')
    #no_dir_mc_early_pulse_mask = dataclasses.I3RecoPulseSeriesMapMask(frame, 'SplitInIcePulses')
    weight_dict = frame['I3MCWeightDict']


    isnuE=False
    isGR=True

    leading_muon_energy = 0

    '''
    reco_particle = frame['EventGeneratorSelectedReco_I3Particle']
    reco_cascade_dir = np.array([reco_particle.dir.x, reco_particle.dir.y, reco_particle.dir.z])
    reco_interaction_pos = reco_particle.pos
    reco_interaction_time = reco_particle.time
    '''
    for particle in mcTree:
        if(particle.type==dataclasses.I3Particle.Hadrons):
            isGR=False
            hadrons=particle
            interaction_time=hadrons.time
            interaction_pos=hadrons.pos
            cascade_dir=np.array([hadrons.dir.x, hadrons.dir.y, hadrons.dir.z])
            break

    for particle in mcTree:
        if particle.type in [dataclasses.I3Particle.ParticleType.MuPlus,
                  dataclasses.I3Particle.ParticleType.MuMinus]:
            muonIds = np.append(muonIds, particle.id.minorID)
            for child in mcTree.children(particle.id):
                muonChildIds = np.append(muonChildIds, child.id.minorID)

    print('muons', muonIds)
    print('muon children', muonChildIds)

    asd

    # If GR, the W boson saved as EMinus in MCTree, so no hadrons
    if isGR==True:
        for particle in mcTree:
            if(particle.type==dataclasses.I3Particle.EMinus):
                hadrons=particle
                interaction_time=hadrons.time
                interaction_pos=hadrons.pos
                cascade_dir=np.array([hadrons.dir.x, hadrons.dir.y, hadrons.dir.z])
                break
    
    for particle in mcTree:
        if(((particle.type==dataclasses.I3Particle.MuMinus) or (particle.type==dataclasses.I3Particle.MuPlus)) and
            (particle.energy>leading_muon_energy)):
            mc_muon = particle
            frame.Put("mc_muon",mc_muon)
            leading_muon_energy=particle.energy
            break
    
    min_distance_to_dom = 2000*I3Units.m

    # RECO
    '''
    reco_dom_distance=np.zeros((87,67))
    reco_dom_hit_times=np.zeros((87,67))
    reco_dom_early_hit_counts=np.zeros((87,67))
    reco_dom_predicted_em_hit_times=np.zeros((87,67))
    reco_dom_predicted_muon_hit_times=np.zeros((87,67))
    reco_dom_hit_sources=np.zeros((87,67))
    reco_cos_separation=np.zeros((87,67))
    reco_max_dom_dist = 0
    '''
    #MC
    dom_distance=np.zeros((87,67))
    dom_hit_times=np.zeros((87,67))
    dom_early_hit_counts=np.zeros((87,67))
    dom_predicted_em_hit_times=np.zeros((87,67))
    dom_predicted_muon_hit_times=np.zeros((87,67))
    dom_hit_sources=np.zeros((87,67))
    cos_separation=np.zeros((87,67))
    max_dom_dist = 0
    #MC_no_dir
    no_dir_max_dom_dist = 0


    for om, geo_info in om_geometry:

        #Removing omkeys from pulse masks
        if(hits.get(om) is not None):
            #reco_early_pulse_mask.set(om, False)
            mc_early_pulse_mask.set(om, False)
            mc_muons_pulse_mask.set(om, False)
            #no_dir_mc_early_pulse_mask.set(om, False)

        #Get all pulses from muons, not necessarily early pulses
        #if (hits.get(om) is not None):

            
        #RECO
        '''
        reco_dot_product = np.dot((geo_info.position-reco_interaction_pos),reco_cascade_dir)
        reco_distance_to_dom = np.sqrt((geo_info.position.x-reco_interaction_pos[0])**2.
                                  +(geo_info.position.y-reco_interaction_pos[1])**2.
                                  +(geo_info.position.z-reco_interaction_pos[2])**2.)
        reco_cos_separation[om[0],om[1]] = reco_dot_product/(reco_distance_to_dom)
        reco_dom_distance[om[0],om[1]]=reco_distance_to_dom
        '''
        #MC
        dot_product = np.dot((geo_info.position-interaction_pos),cascade_dir)
        distance_to_dom = np.sqrt((geo_info.position.x-interaction_pos[0])**2.
                                  +(geo_info.position.y-interaction_pos[1])**2.
                                  +(geo_info.position.z-interaction_pos[2])**2.)
        cos_separation[om[0],om[1]] = dot_product/(distance_to_dom)
        dom_distance[om[0],om[1]]=distance_to_dom

        #Prediction curves
        index_of_ref = 1.33
        #RECO
        '''
        reco_dom_predicted_em_hit_times[om[0],om[1]]=(reco_distance_to_dom/(3e-1/index_of_ref))+reco_interaction_time
        reco_dom_predicted_muon_hit_times[om[0],om[1]]=(reco_distance_to_dom/(3e-1))+reco_interaction_time
        '''
        #MC
        dom_predicted_em_hit_times[om[0],om[1]]=(distance_to_dom/(3e-1/index_of_ref))+interaction_time
        dom_predicted_muon_hit_times[om[0],om[1]]=(distance_to_dom/(3e-1))+interaction_time

        cosThetaC=1/index_of_ref
        sinThetaC=np.sqrt(1-cosThetaC**2)
        #mRange = 800*(1-.5*pow(2.71828,-neutrino_energy/20000)) #correct e
        mRange = 1000
        pRange = 150

        #Saving early hit doms
        #RECO and MC
        #reco_early_dom=np.zeros((0,2))
        mc_early_dom=np.zeros((0,2))
        no_dir_mc_early_dom=np.zeros((0,2))
        
        #updated cone:
        #RECO
        '''
        if (reco_dot_product<(mRange+pRange*cosThetaC) and
            reco_distance_to_dom*np.sqrt(1-reco_cos_separation[om[0],om[1]]**2)<pRange*sinThetaC and 
            reco_cos_separation[om[0],om[1]]>cosThetaC):
            if ((hits.get(om) is not None) and (reco_dom_distance[om[0],om[1]]>reco_max_dom_dist)):
                for reco_hit_counter in range(len(hits.get(om))):
                    if(hits.get(om)[reco_hit_counter].time<reco_dom_predicted_em_hit_times[om[0],om[1]] and 
                       hits.get(om)[reco_hit_counter].time>reco_dom_predicted_muon_hit_times[om[0],om[1]]):
                        reco_max_dom_dist = reco_dom_distance[om[0],om[1]]
                        #Keeping early hits
                        reco_early_pulse_mask.set(om, True)
                        reco_early_dom = np.vstack([reco_early_dom,[om[0],om[1]]])
                        break
        '''
        #MC
        if (dot_product<(mRange+pRange*cosThetaC) and
            distance_to_dom*np.sqrt(1-cos_separation[om[0],om[1]]**2)<pRange*sinThetaC and 
            cos_separation[om[0],om[1]]>cosThetaC):
            if ((hits.get(om) is not None) and (dom_distance[om[0],om[1]]>max_dom_dist)):
                for hit_counter in range(len(hits.get(om))):
                    if(hits.get(om)[hit_counter].time<dom_predicted_em_hit_times[om[0],om[1]] and 
                       hits.get(om)[hit_counter].time>dom_predicted_muon_hit_times[om[0],om[1]]):
                        max_dom_dist = dom_distance[om[0],om[1]]
                        #Keeping early hits
                        mc_early_pulse_mask.set(om, True)
                        mc_early_dom = np.vstack([mc_early_dom,[om[0],om[1]]])
                        break
        '''
        #MC no direction (only condition is distance and 90 degree angular error)
        if ((dom_distance[om[0],om[1]]<750) and cos_separation[om[0],om[1]]>0):
            if ((hits.get(om) is not None) and (dom_distance[om[0],om[1]]>no_dir_max_dom_dist)):
                for hit_counter in range(len(hits.get(om))):
                    if(hits.get(om)[hit_counter].time<dom_predicted_em_hit_times[om[0],om[1]] and 
                       hits.get(om)[hit_counter].time>dom_predicted_muon_hit_times[om[0],om[1]]):
                        no_dir_max_dom_dist = dom_distance[om[0],om[1]]
                        #Keeping early hits
                        print('em time: ', dom_predicted_em_hit_times[om[0],om[1]])
                        print('muon time: ',dom_predicted_muon_hit_times[om[0],om[1]])
                        no_dir_mc_early_pulse_mask.set(om, True)
                        no_dir_mc_early_dom = np.vstack([no_dir_mc_early_dom,[om[0],om[1]]])
                        break
        '''
    
    #RECO
    '''
    if(reco_max_dom_dist>=energy_range[0,1]):
        reco_lmrEnergy = find_nearest(energy_range,reco_max_dom_dist)
    else:
        reco_lmrEnergy = 0
    '''
    #MC
    if(max_dom_dist>=energy_range[0,1]):
        mc_lmrEnergy = find_nearest(energy_range,max_dom_dist)
    else:
        mc_lmrEnergy = 0
    #MC no direction
    if(no_dir_max_dom_dist>=energy_range[0,1]):
        no_dir_mc_lmrEnergy = find_nearest(energy_range,no_dir_max_dom_dist)
    else:
        no_dir_mc_lmrEnergy = 0



    #frame['reco_early_SplitInIcePulses'] = reco_early_pulse_mask
    frame['mc_early_SplitInIcePulses'] = mc_early_pulse_mask
    frame['no_dir_mc_early_SplitInIcePulses'] = no_dir_mc_early_pulse_mask

    #Add muons
    planted_muon = dataclasses.I3Particle()
    planted_muon.speed = dataclasses.I3Constants.c
    planted_muon.location_type = dataclasses.I3Particle.InIce
    planted_muon.shape = dataclasses.I3Particle.Null
    planted_muon.type = dataclasses.I3Particle.MuPlus
    #Reco
    '''
    reco_muon = copy.copy(planted_muon)
    reco_muon.length = reco_max_dom_dist
    reco_muon.pos=reco_interaction_pos
    reco_muon.dir=reco_particle.dir
    reco_muon.time = reco_interaction_time
    reco_muon.energy = reco_lmrEnergy

    reco_early_muon_mcTree.insert(reco_muon)
    frame.Put("reco_early_muon",reco_muon)
    frame.Put("reco_early_muon_mcTree",reco_early_muon_mcTree)
    '''
    #MC (Reco)
    mc_reco_muon = copy.copy(planted_muon)
    mc_reco_muon.length = max_dom_dist
    mc_reco_muon.pos=interaction_pos
    mc_reco_muon.dir=hadrons.dir
    mc_reco_muon.time = interaction_time
    mc_reco_muon.energy = mc_lmrEnergy

    mc_reco_early_muon_mcTree.insert(mc_reco_muon)
    frame.Put("mc_reco_early_muon",mc_reco_muon)
    frame.Put("mc_reco_early_muon_mcTree",mc_reco_early_muon_mcTree)
    #MC no dir (Reco)
    no_dir_mc_reco_muon = copy.copy(planted_muon)
    no_dir_mc_reco_muon.length = no_dir_max_dom_dist
    no_dir_mc_reco_muon.pos=interaction_pos
    no_dir_mc_reco_muon.dir=hadrons.dir
    no_dir_mc_reco_muon.time = interaction_time
    no_dir_mc_reco_muon.energy = no_dir_mc_lmrEnergy

    no_dir_mc_reco_early_muon_mcTree.insert(no_dir_mc_reco_muon)
    frame.Put("no_dir_mc_reco_early_muon",no_dir_mc_reco_muon)
    frame.Put("no_dir_mc_reco_early_muon_mcTree",no_dir_mc_reco_early_muon_mcTree)

    '''
    with open("/data/user/eyildizci/leadingMuonReco/lmr_"+str(mRange)+"_"+str(pRange)+"_"+str(infileNumber)+".txt", "a") as myfile:
        myfile.write(str(reco_lmrEnergy) + '   ' + str(lmrEnergy) + '   ' + str(leading_muon_energy) + '\n')
    '''

#Run the I3Tray object, begins loop over frames
tray.AddModule(EarlyHitIdentifier,'earlyhitidentifier',Streams = [icetray.I3Frame.Physics])
tray.AddModule("I3Writer", Filename='./i3files/early_hits_PEPE_'+str(infileNumber)+'.i3.zst',
               DropOrphanStreams=[icetray.I3Frame.TrayInfo],
               Streams=[icetray.I3Frame.TrayInfo,
                        icetray.I3Frame.DAQ,
                        icetray.I3Frame.Physics])   
tray.Execute() #can specify number of frames you want to run on
#closes I3Tray
tray.Finish()