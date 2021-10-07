import ROOT, numpy as np, pandas as pd, sys, root_pandas,time
ROOT.gSystem.Load('../delphes_install/lib/libDelphes')

n_bad_tracks = 0
n_tracks = 0

def JB(event):

    VAP = 0
    VP = 0
    branchEFlowTrack = event.EFlowTrack
    branchElectron = event.Electron
    branchEFlowPhoton = event.EFlowPhoton
    branchEFlowNeutralHadron = event.EFlowNeutralHadron
    delta_track = 0
    temp_p = ROOT.TVector3(0,0,0)
    for i in range(branchEFlowTrack.GetEntries()):
       track_mom = branchEFlowTrack.At(i).P4()
       #print(track_mom.Pz(),track_mom.E())
       delta_track += (track_mom.E() - track_mom.Pz())
       temp_p = temp_p + track_mom.Vect()
       
       #debug: check fraction of bad tracks
       global n_bad_tracks
       global n_tracks
       if np.isnan(track_mom.Pz()):
           n_bad_tracks+=1
       n_tracks+=1
       print("bad track fraction",n_bad_tracks/n_tracks)
       #print(track_mom.E() , ' ' , track_mom.Pz())
    if branchElectron.GetEntries()>0:
        e = branchElectron.At(0).P4()
        delta_track_noel = delta_track - (e.E() - e.Pz())
    delta_photon = 0
    for i in range(branchEFlowPhoton.GetEntries()):
       pf_mom = branchEFlowPhoton.At(i).P4()
       delta_photon += (pf_mom.E() - pf_mom.Pz())
       temp_p = temp_p + pf_mom.Vect()
    delta_neutral = 0
    delta_neutral_noBarrel = 0
    pt_noBarrel = ROOT.TVector3()
    for i in range(branchEFlowNeutralHadron.GetEntries()):
       pf_mom = branchEFlowNeutralHadron.At(i).P4()
       delta_neutral += (pf_mom.E() - pf_mom.Pz())
       temp_p = temp_p+ pf_mom.Vect()
       if abs(pf_mom.Eta())>1.0:
           delta_neutral_noBarrel += (pf_mom.E() - pf_mom.Pz())
           pt_noBarrel = pt_noBarrel + pf_mom.Vect()
    delta = delta_track + delta_photon + delta_neutral
    
    
    
    
    
    #if np.isnan(delta) :
    #print(delta, delta_track, delta_photon, delta_neutral)
    
    y_JB   = (delta)/(2.0*10.0)
    ptmiss = temp_p.Perp()
    
    for i in range(branchEFlowPhoton.GetEntries()):
        pf_mom = branchEFlowPhoton.At(i).P4()
        dot =pf_mom.X() * temp_p.X() + pf_mom.Y() * temp_p.Y()
        if(dot > 0):
            VP += dot
        else :
            VAP += -dot
            
        
    for i in range(branchEFlowNeutralHadron.GetEntries()):
        pf_mom = branchEFlowNeutralHadron.At(i).P4()
        dot =pf_mom.X() * temp_p.X() + pf_mom.Y() * temp_p.Y()
        if(dot > 0):
            VP += dot
        else :
            VAP += -dot
    
    
    Q2_JB  = (ptmiss*ptmiss)/(1-y_JB)
    
    s     = 4*10.0*275.0
    if(y_JB>0):
        x_JB  = Q2_JB/(s*y_JB)
    else:
        x_JB = -999
    pt_all = temp_p
    return pt_all.Pt(), pt_all.Eta(), pt_all.Phi(), Q2_JB, x_JB,y_JB,delta,VAP, VP

def convert(infilename, outfilename,debug=False,N=None,hadronTuple=False,arg_maxR = 0.9):
    start = time.perf_counter()
    infile = ROOT.TFile(infilename)
    tree = infile.Delphes
    d = {}
    #fields = "PT Eta Phi Mass DeltaEta DeltaPhi Flavor FlavorAlgo FlavorPhys BTag BTagAlgo BTagPhys TauTag TauWeight Charge EhadOverEem NCharged NNeutrals NeutralEnergyFraction ChargedEnergyFraction NSubJetsTrimmed NSubJetsPruned NSubJetsSoftDropped".split()
    fields = "PT Eta Phi Mass DeltaEta DeltaPhi NCharged NNeutrals NeutralEnergyFraction ChargedEnergyFraction".split()
    jetFields = ["Jet_" + name for name in fields]
    genJetFields = ["GenJet_" + name for name in fields]
    
    
    #fields = "PID Status IsPU M1 M2 D1 D2 Mass E Px Py Pz PT Eta Phi Rapidity D0 DZ".split()
    fields = "PID Status E Px Py Pz PT Eta".split()
    
    neutrinoFields = ["Neutrino_" + f for f in fields]
    quarkFields = ["Quark_" + f for f in fields]
    
    
    fields = "PID Charge P PT Eta Phi E ET Eem Ehad".split()
    hadronFields = ["Hadron_" + f for f in fields]
    fields = "PID Charge P PT Eta Phi E".split()
    hadronFields += ["GenHadron_" + f for f in fields]
    
    otherFields = ["MissingET_MET", "MissingET_Eta", "MissingET_Phi","GenMissingET_MET", "GenMissingET_Eta","GenMissingET_Phi","Event_Number"]

    allFields = jetFields + genJetFields + neutrinoFields + quarkFields + otherFields + "Gen_W2 Gen_x Gen_y Gen_Q2".split()

    if hadronTuple:
        allFields+= hadronFields
        # Index of hadron within event.
        # To create a subtuple with exactly one entry per event, one can do df.query("Hadron_i == 0")
        allFields.append("Hadron_i")

    # index of the jet within the event.
    # To create a subtuple with exactly one entry per event, one can do df.query("Jet_i == 0")
    allFields.append("Jet_i")
    
    #set up counters for every type of common constituent (>1% of the total constituents.
    pids = (22,211,-211,2212,-2212,2112,-2112,130,310,321,-321)
    for pid in pids :
        allFields.append("GenJet_n_"+("m" if pid <0 else "") + str(abs(pid)))
    allFields.append("GenJet_n_other")
    
    #modify the names of the columns in the output (so that it's easier to use)
    inputNames = {name : name.replace("_",".").replace("Neutrino","Particle").replace("Quark","Particle").replace("Hadron_","") for name in allFields}
    for name in allFields:
        d[name] = []


    for i,event in enumerate(tree):
        if N != None and i>N and N>0:
            break
        
        match_indices = []
        
        met_JB,eta_JB,phi_JB, Q2_JB, x_JB,y_JB,delta_JB,VAP,VP = JB(event)
        if not "JB_MET" in d.keys():
            for f in 'JB_MET JB_Eta JB_Phi JB_Q2 JB_x JB_y JB_delta VAP VP'.split():
                d[f] = []
        
        
        nu_properties = {}
        q_properties = {}
        other_properties = {}
        #find neutrino and quark
        
        
        for particle in event.Particle:
            pid = particle.PID
            status = particle.Status
            if pid == 12 and status == 1:
                for name in neutrinoFields:
                    #print(name)
                    nu_properties[name] = getattr(particle,name.replace("Neutrino_",""))
            elif status == 23 and pid in [1,3,5, -2, -4]:
                for name in quarkFields:
                    #print(name)
                    q_properties[name] =  getattr(particle,name.replace("Quark_",""))
        
        branchParticle=event.Particle
        pProton      = branchParticle.At(0).P4(); #these numbers 0 , 3, 5 are hardcoded in Pythia8
        pleptonIn    = branchParticle.At(3).P4();
        pleptonOut   = branchParticle.At(5).P4();
        pPhoton      = pleptonIn - pleptonOut;
        #Q2, W2, Bjorken x, y, nu.
        Q2 = -pPhoton.M2()
        W2 = (pProton + pPhoton).M2()
        x = Q2 / (2. * pProton.Dot(pPhoton))
        y = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn))
        other_properties['Gen_Q2'] = Q2
        other_properties['Gen_W2'] = W2
        other_properties['Gen_x'] = x
        other_properties['Gen_y'] = y
        
        njets = int(tree.GetLeaf("Jet_size").GetValue())
        ngenjets = int(tree.GetLeaf("GenJet_size").GetValue())
        
        
        for name in otherFields:
            other_properties[name] = tree.GetLeaf(inputNames[name]).GetValue(0)
        if i%1000 == 0 and i != 0:
            print(other_properties["Event_Number"], "events,", (time.perf_counter()-start), " s,   (",(time.perf_counter()-start)/other_properties["Event_Number"], " s avg)")
        #first find matches
        for jet in event.Jet:
            eta = jet.Eta
            phi = jet.Phi
            maxR = arg_maxR
            kbest = -1
            Rbest = 9999999
            for k,genJet in enumerate(event.GenJet):
                etagen = genJet.Eta
                phigen = genJet.Phi
                dphi = phi-phigen
                if dphi > np.pi:
                    dphi -= phi
                if dphi < -np.pi:
                    dphi += phi
                R = np.hypot(eta-etagen,dphi)
                if R< Rbest and R<maxR:
                    Rbest = R
                    kbest = k
            match_indices.append(kbest)
        if(debug):
            print("finding matches complete")
        
        #count the number of each type of constituent in the generated jets:
        counts = {}
        for pid in pids:
            counts[pid] = []
        count_other = []
        for genjet in event.GenJet:
            count_other.append(0)
            for pid in pids:
                counts[pid].append(0)
            for particle in genjet.Particles:
                #print("constituent pid = " + str(particle.PID))
                if int(particle.PID) in pids:
                    pid = particle.PID
                    counts[pid][-1] = 1+counts[pid][-1]
                    
                else:
                    count_other[-1] = 1+count_other[-1]
        
        #loop through recon jets
        for j,jet in enumerate(event.Jet):
            ntracks = 0
            for h,(track,particle) in enumerate(zip(jet.Constituents,jet.Particles)):
                #if not "Track" in str(type(track)):
                #    continue
                for name in jetFields:
                    d[name].append(getattr(jet,name.replace("Jet_","")))
                for name in genJetFields:
                    if(debug):
                        print(name)
                    if match_indices[j] != -1:
                        d[name].append(tree.GetLeaf(inputNames[name]).GetValue(match_indices[j]))
                    else :
                        d[name].append(0)
                if match_indices[j] != -1:
                    for pid in pids:
                        fieldname = "GenJet_n_"+("m" if pid <0 else "") + str(abs(pid))
                        d[fieldname].append(counts[pid][match_indices[j]])
                    d["GenJet_n_other"].append(count_other[match_indices[j]])
                else :
                    for pid in pids:
                        fieldname = "GenJet_n_"+("m" if pid <0 else "") + str(abs(pid))
                        d[fieldname].append(0)
                    d["GenJet_n_other"].append(0)
                
                for name in neutrinoFields:
                    d[name].append(nu_properties[name])
                for name in quarkFields:
                    d[name].append(q_properties[name])
                for name in other_properties.keys():
                    d[name].append(other_properties[name])
                d['Jet_i'].append(j)
                
                
                #there are some fields that only apply to tracks, and others that only apply to towers.
                # 0 represents N/A
                if "Track" in str(type(track)):
                    for name in hadronFields:
                        if not 'Gen' in name:
                            if name not in "Hadron_E Hadron_ET Hadron_Eem Hadron_Ehad".split():
                                d[name].append(getattr(track,name.replace("Hadron_","")))
                            else:
                                d[name].append(0)
                        else:
                            d[name].append(getattr(particle,name.replace("GenHadron_","")))
                else :
                    for name in hadronFields:
                        if not 'Gen' in name:
                            if name not in "Hadron_P Hadron_PT Hadron_PID Hadron_Charge".split():
                                d[name].append(getattr(track,name.replace("Hadron_","")))
                            else:
                                d[name].append(0)
                        else:
                            d[name].append(getattr(particle,name.replace("GenHadron_","")))
                d['Hadron_i'].append(ntracks)
                d['JB_MET'].append(met_JB)
                d['JB_Eta'].append(eta_JB)
                d['JB_Phi'].append(phi_JB)
                d['JB_Q2'].append(Q2_JB)
                d['JB_x'].append(x_JB)
                d['JB_y'].append(y_JB)
                d['JB_delta'].append(delta_JB)
                d['VAP'].append(VAP)
                d['VP'].append(VP)
                ntracks+=1
            # jet without constituents (does this ever happen?)
            if ntracks == 0:
                for name in jetFields:
                    d[name].append(getattr(jet,name.replace("Jet_","")))
                for name in genJetFields:
                    if(debug):
                        print(name)
                    if match_indices[j] != -1:
                        d[name].append(tree.GetLeaf(inputNames[name]).GetValue(match_indices[j]))
                    else :
                        d[name].append(0)
                for name in neutrinoFields:
                    d[name].append(nu_properties[name])
                for name in quarkFields:
                    d[name].append(q_properties[name])
                for name in other_properties.keys():
                    d[name].append(other_properties[name])
                d['Jet_i'].append(j)
                for name in hadronFields:
                    d[name].append(0)
                
                d['Hadron_i'].append(0)
                d['JB_MET'].append(met_JB)
                d['JB_Eta'].append(eta_JB)
                d['JB_Phi'].append(phi_JB)
                d['JB_Q2'].append(Q2_JB)
                d['JB_x'].append(x_JB)
                d['JB_y'].append(y_JB)
                d['JB_delta'].append(delta_JB)
                d['VAP'].append(VAP)
                d['VP'].append(VP)
                if match_indices[j] != -1:
                    for pid in pids:
                        fieldname = "GenJet_n_"+("m" if pid <0 else "") + str(abs(pid))
                        d[fieldname].append(counts[match_indices[j]][pid])
                    d["GenJet_n_other"].append(count_other[match_indices[j]])
                else :
                    for pid in pids:
                        fieldname = "GenJet_n_"+("m" if pid <0 else "") + str(abs(pid))
                        d[fieldname].append(0)
                    d["GenJet_n_other"].append(0)
        if(debug):
            print("adding recon jets complete")
        
        #now add unmatched gen jets:
        ii = njets
        for k,genJet in enumerate(event.GenJet):
            if not k in match_indices:
                for name in genJetFields:
                    d[name].append(getattr(genJet,name.replace("GenJet_","")))
                for name in jetFields:
                    d[name].append(0)
                for name in neutrinoFields:
                    d[name].append(nu_properties[name])
                for name in quarkFields:
                    d[name].append(q_properties[name])
                for name in other_properties.keys():
                    d[name].append(other_properties[name])
                for name in hadronFields:
                    d[name].append(0)
                d['Hadron_i'].append(0)
                d['Jet_i'].append(ii)
                d['JB_MET'].append(met_JB)
                d['JB_Eta'].append(eta_JB)
                d['JB_Phi'].append(phi_JB)
                d['JB_Q2'].append(Q2_JB)
                d['JB_x'].append(x_JB)
                d['JB_y'].append(y_JB)
                d['JB_delta'].append(delta_JB)
                d['VAP'].append(VAP)
                d['VP'].append(VP)
                for pid in pids:
                    fieldname = "GenJet_n_"+("m" if pid <0 else "") + str(abs(pid))
                    d[fieldname].append(counts[pid][k])
                d["GenJet_n_other"].append(count_other[k])
                ii+=1
        if(debug):
            print("adding unmatched gen jets complete")
        #if there are no gen jets nor recon jets, add an entry for the event itself
        if ngenjets == 0 and njets == 0:
            for name in genJetFields:
                d[name].append(0)
            for name in jetFields:
                d[name].append(0)
            for name in neutrinoFields:
                d[name].append(nu_properties[name])
            for name in quarkFields:
                d[name].append(q_properties[name])
            for name in other_properties.keys():
                d[name].append(other_properties[name])
            for name in hadronFields:
                d[name].append(0)
            d['Hadron_i'].append(0)
            d['Jet_i'].append(0)
            d['JB_MET'].append(met_JB)
            d['JB_Eta'].append(eta_JB)
            d['JB_Phi'].append(phi_JB)
            d['JB_Q2'].append(Q2_JB)
            d['JB_x'].append(x_JB)
            d['JB_y'].append(y_JB)
            d['JB_delta'].append(delta_JB)
            d['VAP'].append(VAP)
            d['VP'].append(VP)
            for pid in pids:
                fieldname = "GenJet_n_"+("m" if pid <0 else "") + str(abs(pid))
                d[fieldname].append(0)
            d["GenJet_n_other"].append(0)
        del q_properties
        del nu_properties
        
        
        
        
    print({a: len(d[a]) for a in d.keys() })
    df = pd.DataFrame(d)
    
    #save some space by converting int64 fields to int16 fields:
    ints = df.select_dtypes(include=['int64']).columns.tolist()
    df[ints] = df[ints].astype('int16')
    
    '''
    df['Neutrino_Pz_fromMiss'] = df.eval("MissingET_MET*sinh(MissingET_Eta)+265")
    df['Neutrino_E_fromMiss']  = df.eval("sqrt((MissingET_MET*sinh(MissingET_Eta)+265)**2 + MissingET_MET**2)")
    df['Neutrino_Pz_fromGenMiss'] = df.eval("GenMissingET_MET*sinh(GenMissingET_Eta)+265")
    df['Neutrino_E_fromGenMiss']  = df.eval("sqrt((GenMissingET_MET*sinh(GenMissingET_Eta)+265)**2 + GenMissingET_MET**2)")
    df['Q2'] = df.eval('2*10*(Neutrino_E_fromMiss+Neutrino_Pz_fromMiss)')
    df['Xb'] = df.eval('Q2/(2*((10-Neutrino_E_fromMiss)*275-(-10-Neutrino_Pz_fromMiss)*sqrt(275**2-.9383**2)))')
    df['ye'] = df.eval('((10-Neutrino_E_fromMiss)*275-(-10-Neutrino_Pz_fromMiss)*sqrt(275**2-.9383**2))/(2*10*275)')
    '''
    
    df["Jet_E"] = df.eval("sqrt((Jet_PT*cosh(Jet_Eta))**2+Jet_Mass**2)")
    df["GenJet_E"] = df.eval("sqrt((GenJet_PT*cosh(GenJet_Eta))**2+GenJet_Mass**2)")
    #df["QuarkI_E"] = df.eval("Quark_E+(-10+Neutrino_E)")
    #df["QuarkI_Pz"] = df.eval("Quark_Pz+(10+Neutrino_Pz)")
    #df["QuarkI_Px"] = df.eval("Quark_Px+(Neutrino_Px)")
    #df["QuarkI_Py"] = df.eval("Quark_Py+(Neutrino_Py)")
    #df["QuarkI_Phi"] = np.arctan2(df.QuarkI_Py,df.QuarkI_Px)
    #df["QuarkI_PT"] = df.eval("sqrt(QuarkI_Py**2+QuarkI_Px**2)")


    df['Hadron_Px'] = df.eval('Hadron_PT*cos(Hadron_Phi)')
    df['Hadron_Py'] = df.eval('Hadron_PT*sin(Hadron_Phi)')
    df['Hadron_Pz'] = df.eval('Hadron_PT*sinh(Hadron_Eta)')
    df['Hadron_P'] = df.eval('Hadron_PT*cosh(Hadron_Eta)')
    
    df['GenHadron_Px'] = df.eval('GenHadron_PT*cos(GenHadron_Phi)')
    df['GenHadron_Py'] = df.eval('GenHadron_PT*sin(GenHadron_Phi)')
    df['GenHadron_Pz'] = df.eval('GenHadron_PT*sinh(GenHadron_Eta)')
    df['GenHadron_P'] = df.eval('GenHadron_PT*cosh(GenHadron_Eta)')

    df['Jet_Px'] = df.eval('Jet_PT*cos(Jet_Phi)')
    df['Jet_Py'] = df.eval('Jet_PT*sin(Jet_Phi)')
    df['Jet_Pz'] = df.eval('Jet_PT*sinh(Jet_Eta)')
    df['Jet_P'] = df.eval('Jet_PT*cosh(Jet_Eta)')

    df['GenJet_Px'] = df.eval('GenJet_PT*cos(GenJet_Phi)')
    df['GenJet_Py'] = df.eval('GenJet_PT*sin(GenJet_Phi)')
    df['GenJet_Pz'] = df.eval('GenJet_PT*sinh(GenJet_Eta)')
    df['GenJet_P'] = df.eval('GenJet_PT*cosh(GenJet_Eta)')

    expr = "sqrt(GenMissingET_MET**2+GenJet_PT**2+2*GenMissingET_MET*GenJet_PT*cos(GenJet_Phi-GenMissingET_Phi))"
    df['Gen_qT'] = df.eval(expr)
    df['qT'] = df.eval(expr.replace("Gen",""))

    expr = "arctan2(GenMissingET_MET*sin(GenMissingET_Phi)+GenJet_Py,GenMissingET_MET*cos(GenMissingET_Phi)+GenJet_Px)"
    df['Gen_qT_Phi'] = df.eval(expr)
    df['qT_Phi'] = df.eval(expr.replace("Gen",""))


    df['Hadron_Jt'] = df.eval("sqrt((Hadron_Px*Jet_Py-Hadron_Py*Jet_Px)**2+"\
                            + "(Hadron_Pz*Jet_Px-Hadron_Px*Jet_Pz)**2+"\
                            +  "(Hadron_Px*Jet_Py-Hadron_Py*Jet_Px)**2)/Jet_P")
    df['Hadron_Zh'] = df.eval('(Hadron_Px*Jet_Px+Hadron_Py*Jet_Py+Hadron_Pz*Jet_Pz)/Jet_P**2')
    
    df['GenHadron_Jt'] = df.eval("sqrt((GenHadron_Px*GenJet_Py-GenHadron_Py*GenJet_Px)**2+"\
                            + "(GenHadron_Pz*GenJet_Px-GenHadron_Px*GenJet_Pz)**2+"\
                            +  "(GenHadron_Px*GenJet_Py-GenHadron_Py*GenJet_Px)**2)/GenJet_P")
    df['GenHadron_Zh'] = df.eval('(GenHadron_Px*GenJet_Px+GenHadron_Py*GenJet_Py+GenHadron_Pz*GenJet_Pz)/GenJet_P**2')
    
    df.to_root(outfilename, "jets")

if __name__ == '__main__':
    infilename=sys.argv[1]
    outfilename=sys.argv[2]
    n = None
    if len(sys.argv)>2:
        n = int(sys.argv[3])
    if "-h" in sys.argv:
        hadronTuple = True
    else:
        hadronTuple = False
    convert(infilename,outfilename, N=n,hadronTuple=hadronTuple,arg_maxR=0.9)

