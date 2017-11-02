#TODO: implement c-jet and l-jet control regions
import os,sys
sys.path.append('/uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ZplusC/python/')
import ConfigParser
import myutils as utils 
import ROOT,array

ROOT.gSystem.Load('/uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ZplusC/interface/VHbbNameSpace_h.so')

ROOT.gROOT.SetBatch(True)

def AddFileToChain(chain, sampleName, cfg, path, debug = False):
    if debug: print '>>>>>>>Sample used: >>>>>>>>>>>>'
    for i in range(int(cfg.get('General','nSamples'))):
      sample = 'Sample_' + str(i)
      name = cfg.get(sample,'name')
      fileName = path + '/' + cfg.get(sample,'file')
      if name == sampleName: #add all dataset with the same name together
        if debug: print fileName
        chain.Add(fileName)
        if debug: print 'Current chain events: ', chain.GetEntries()

def make1Dplot(hDir, h2D, varN, bins): #parent directory, 2D histogram, pt_binning (string)
  #ProjectionX(name,firstbin,lastbin)
  #loop over bins
  #-make the projected plot name
  #-save plot to folder
  for iBin in range(len(bins)-1):
      binName = varN + '_' + bins[iBin] + '_' + bins[iBin+1]
      iBin_h2D = h2D.GetYaxis().FindFixBin((float(bins[iBin]) + float(bins[iBin+1]))/2.)
      hDir.cd(binName)
      plotName = h2D.GetName()
      plotName = plotName.replace('2D',binName.replace(varN + '_',''))
      hP = h2D.ProjectionX(plotName, iBin_h2D, iBin_h2D)
      hP.Write()
      hDir.cd()
'''
def makeVarForDraw(varName, jetIdx= ''):
    if varName == 'Jet_pt':
        if jetIdx == '': raise ValueError('Empty jet index')
        return 'Jet_pt['+jetIdx+']'
    elif varName == 'V_pt':
        return 'VHbb::Vpt(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0],vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1])'
    else: raise ValueError('Undefined variable name. Check configuration')
    return ''
'''
def makeVar(varIn,jetIdxIn):
  if 'idxJet' in varIn: 
    return varIn.replace('idxJet',jetIdxIn)
  return varIn
    

def getVarAttribute(varName, varAtt):
    atts = varAtt.replace(' ','').split('},{')
    return [att.replace('{','').replace('}','').replace(varName + ':','') for att in atts if varName in att][0]
    
def makeBinning(nBin,xL,xH):
  binningTmp = []
  nBin = int(nBin)
  xL = float(xL)
  xH = float(xH)
  bW = (xH-xL)/nBin
  for iB in range(nBin):
    binningTmp.append(xL + iB*bW)
  binningTmp.append(xH)
  return binningTmp


#####################
#Main
#####################

#####################
#some settings
#####################
cfg = utils.BetterConfigParser()
cfg.read('Config/config.ini')

#>>>>>>>>>>>>>>>Modify here>>>>>>>>>>>>
useDYnlo = False
fName = 'Test/test_histo_lepPt_25_pt_eta_3D_plots_binning_v2_xxx.root'
#regions = ['zjet','zHFjet','emu','ttsemi'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
#regions = ['zHFjet','emu','ttsemi'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
regions = ['zHFjet','ttsemi'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
#regions = ['wjet'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
#labels used in dictionary which store TChains of data files

dataLabels = {'zjet':{'Data_Ele':'DoubleEG','Data_Muon':'DoubleMuon','DY':'DY','TT':'TT','WW':'WW','WZ':'WZ','ZZ':'ZZ'}, #region,data for each channel, MC samples 
              'zHFjet':{'Data_Ele':'DoubleEG','Data_Muon':'DoubleMuon','DY':'DY','TT':'TT','WW':'WW','WZ':'WZ','ZZ':'ZZ'},
              'emu':{'Data_Muon':'SingleMuon','DY':'DY','TT':'TT','WW':'WW','WZ':'WZ','ZZ':'ZZ'},
              'ttsemi':{'Data_Muon':'SingleMuon','DY':'DY','TT':'TT','WW':'WW','WZ':'WZ','ZZ':'ZZ'},
              'wjet':{'Data_Muon':'SingleMuon','Wjet':'Wjet','DY':'DY','TT':'TT','WW':'WW','WZ':'WZ','ZZ':'ZZ'}
             } #label to get dataset

if useDYnlo:
    for k,v in dataLabels.items():
        for k1,v1 in v.items():
            if k1 == 'DY': dataLabels[k][k1] = 'DY_nlo'

print dataLabels

if useDYnlo:
    fName = fName.replace('.root','_DY_nlo.root')

fOut = ROOT.TFile.Open(fName,'recreate')
######################
#get data and MC files and sample cat cuts
######################

#TEMP
#sys.exit()
var_for_2D = cfg.get('General','var_for_2D').split(',')
var_for_2D_binning = cfg.get('General','var_for_2D_binning')
var_for_2D_formula = cfg.get('General','var_for_2D_formula')
varBinnings = {}
varFormulas = {}
for var in var_for_2D:
    val1 = getVarAttribute(var,var_for_2D_binning).split(',') 
    val2 = array.array('d',[])
    for i in val1:
        val2.append(float(i))
    nBins_tmp = len(val2) - 1
    varBinnings[var] = [int(nBins_tmp),val1,val2] #val1 = string, val2 = float
    varFormulas[var] = getVarAttribute(var,var_for_2D_formula)

var_for_3D = cfg.get('General','var_for_3D').split(',')
var_for_3D_binning = cfg.get('General','var_for_3D_binning')
var_for_3D_formula = cfg.get('General','var_for_3D_formula')
var_3D_binnings = {}
var_3D_formulas = {}
for var in var_for_3D:
    val1 = getVarAttribute(var,var_for_3D_binning).split(',') 
    val2 = array.array('d',[])
    for i in val1:
        val2.append(float(i))
    nBins_tmp = len(val2) - 1
    var_3D_binnings[var] = [int(nBins_tmp),val1,val2] #val1 = string, val2 = float
    var_3D_formulas[var] = getVarAttribute(var,var_for_3D_formula)

print '>>>>>>>>>>>>>>>>>>>>'
print var_for_3D
print var_3D_binnings
print var_3D_formulas


jet_vtxMass_name = cfg.get('General','jet_vtxMass_name')

#loop over regions
for r in regions:
  print '##################################'
  print 'Region: ', r
  print '##################################'
  hDir = fOut.mkdir(r)
  hDir.cd()
  #get chain
  path = cfg.get('Paths','path_' + r)
  chains = {}
  for dataName,sampleName in dataLabels[r].items():
    chains[dataName] = ROOT.TChain('tree')
    AddFileToChain(chains[dataName],sampleName,cfg,path)

  #Print number of entries
  print '# Dataset summary: '
  print 'Path: ', path
  for k,v in chains.items():
    print k, ' ', v.GetEntries()
  
  #loop over channel
  for chan in ['Ele','Muon']:

#>>>>>>>>>> MODIFY (CHANGE) HERE>>>>>>>>>>>>>

    if (r == 'emu' or r == 'ttsemi' or r == 'wjet') and chan == 'Ele': continue
    
    rName = r + '_' + chan
    plots = cfg.get('Plots',r + '_plot').split(',')
    cutCommon = cfg.get('Cuts', rName) #get cut for a region
    jsonCut = cfg.get('Cuts','jsonCut')
    sf = cfg.get('SFs','sf_'+rName)
    cat_cuts = cfg.get('Cuts',r+'_cat').split(',')
    cat_cuts.append('') #empty == all jets
    print '\n#######################################'
    print '####Region, chan: ' + r + ' ' + chan
    print 'Common setting'
    print 'Common cuts:  ', cutCommon
    print 'Json cuts:    ', jsonCut
    print 'SFs:   ', sf
    print 'Plots: ', plots
    
    #make directory
    hDir1 = hDir.mkdir(chan)
    #change directory
    hDir1.cd()
    
    #make directory for binned plots
    for varN,varB in varBinnings.items():

#>>>>>>>MODIFY (CHANGE HERE)>>>>>>>>>>>
#only make folder for jet_pt when the region is NOT zHFjet or zjet
        if (r != 'zjet' and r != 'zHFjet') and varN != var_for_2D[0]: continue
        for iBin in range(len(varB[1]) - 1):
            binName = varN + '_' + varB[1][iBin] + '_' + varB[1][iBin+1] #0 = nBins
            hDir1.mkdir(binName)
    
    
    #Save cut for this region and channel
    setup_str = ROOT.TObjString('Common_cuts:' + cutCommon)
    setup_str.Write()
    setup_str = ROOT.TObjString('Json cuts:' + jsonCut)
    setup_str.Write()
    setup_str = ROOT.TObjString('Scale:' + sf)
    setup_str.Write()
    strTmp = 'Jet_cat_cut:' 
    for tmp in cat_cuts:
        strTmp = strTmp + ',' + tmp
    setup_str = ROOT.TObjString(strTmp)
    setup_str.Write()
    
    #loop over samples
    setup_final_str = '' 
    for k,v in chains.items():
      #Temp
      #if 'DY' not in k: continue
      if 'Data' in k and chan not in k: continue #for data only loop over corresponding data set 
      #add json requirement to cut
      cut = cutCommon
      if 'Data' in k: cut = jsonCut + '&&' + cut
      #if 'Data' not in k: cut = '(' + cut + ')' + '*(' + sf + ')'
      print '>>>>>>>>>>>>>>>>>>>>>>>'
      print 'Chan, region, sample: ', chan, ' ', r, ' ', k
      print 'Entries: ', v.GetEntries()
      print 'Cut: ', cut
      for plot in plots:
        
        print '>>>>', plot
        
#>>>>>>>>>>>>>MODIFY (CHANGE) HERE>>>>>>>>>>>>>>

        #plot variable
        jetIdxTmp = ''
        if r != 'ttsemi': jetIdxTmp = cfg.get('General','idxJet_'+r)
        if r == 'ttsemi': jetIdxTmp = cfg.get('General','idxJet_2_ttsemi') 
        
        var = cfg.get(plot,'var')

        #define jet index
        var = makeVar(var,jetIdxTmp)
        print 'Plot variable after define idxJet: ', var
        
        #define x-axis range
        xAxisRange = cfg.get(plot,'range').split(',')
        xAxisRange1 = [] #make conversion to accommodate special case of Pi()
        for i in range(len(xAxisRange)):
          if 'TMath::Pi' in xAxisRange[i] and '-' in xAxisRange[i]:
            xAxisRange1.append(-ROOT.TMath.Pi())
          elif 'TMath::Pi' in xAxisRange[i] and '-' not in xAxisRange[i]:
            xAxisRange1.append(ROOT.TMath.Pi())
          else:
            xAxisRange1.append(xAxisRange[i])
        
        xAxisRange = cfg.get(var_for_2D[0],'range').split(',')
        xAxisRange2 = [xAxisRange[i] for i in range(len(xAxisRange))] #for making finer bin jet_pt plot
        

        sf_final = '(1)'
        if 'Data' not in k: sf_final = sf

#>>>>>>>>>>>>>>MODIFY (CHANGE) here>>>>>>>>>>>>>>>>      
        for cat_cut in cat_cuts:
          
          if 'Data' in k and cat_cut != '': continue #skip making b-,c-,l-jet in case of data
          
          #if 'DY' not in k and (r == 'zjet' or r == 'zHFjet') and cat_cut != '': continue
          if 'DY' not in k and (r == 'zjet') and cat_cut != '': continue
          if 'TT' not in k and (r == 'emu' or r == 'ttsemi') and cat_cut != '': continue
          if 'Wjet' not in k and (r == 'wjet') and cat_cut != '': continue
          
          #make cut
          cutTmp = cut
          if cat_cut != '':
            cutTmp = cut+'&&('+cat_cut.split(':')[1]+')'
          
          #make 1D plot
          plotName = plot + '_' + k.replace('_' + chan,'') + '_' + rName
          if cat_cut != '':
            plotName = plot + '_' + k.replace('_' + chan,'') + '_' + cat_cut.split(':')[0] + '_' + rName
          axis = {'x':[int(xAxisRange1[0]),float(xAxisRange1[1]),float(xAxisRange1[2])]}
          h = utils.util_funcs.getHisto(v, str(plotName), axis, [var], cutTmp, sf_final)
          h.Write()
    

          #make all plots in parametrization bins for zjet and zHFjet region and only make plots for jet_pt parametrization for other regions
          for para_var_name,para_var_binning in varBinnings.items():
            if r == 'zjet' or r == 'zHFjet' or para_var_name == var_for_2D[0]:
                plotName = plot + '_' + para_var_name + '_2D_' + k.replace('_' + chan,'') + '_' + rName
                if cat_cut != '':
                  plotName = plot + '_' + para_var_name + '_2D_' + k.replace('_' + chan,'') + '_' + cat_cut.split(':')[0] + '_' + rName
                  
                axis = {'x':[int(xAxisRange1[0]), float(xAxisRange1[1]), float(xAxisRange1[2])],'y':[para_var_binning[0],array.array('d',para_var_binning[2])]}
                var_for_draw = makeVar(varFormulas[para_var_name],jetIdxTmp) 
                h_2D = utils.util_funcs.getHisto(v, plotName, axis, [var,var_for_draw], cutTmp, sf_final)
                make1Dplot(hDir1,h_2D,para_var_name,para_var_binning[1])
                h_2D.Write()
            
          #make 3D plot of x=jet_vtxMass y=jet_pt and z=jet_eta:
          #1. zHFjet: make this plot for each parameterized variables for example jet_pt V_pt (only make this for MC samples)
          #2. for emu and other validation region: make this plot for all jets disregard of their pt and V_pt (make this for data as well)
          #only do this for jet_vtxMass and cat_cut != '' to do jet_pt and jet_eta reweighting
          if plot == jet_vtxMass_name:
            plotName = plot + '_' + var_for_3D[0] + '_' + var_for_3D[1] + '_3D_' + k.replace('_' + chan,'') + '_' + cat_cut.split(':')[0] + '_' + rName
            if cat_cut == '':
              plotName = plot + '_' + var_for_3D[0] + '_' + var_for_3D[1] + '_3D_' + k.replace('_' + chan,'') + '_' + rName
            axis = {'x':[int(xAxisRange1[0]), float(xAxisRange1[1]), float(xAxisRange1[2])],'y':[int(var_3D_binnings[var_for_3D[0]][2][0]),float(var_3D_binnings[var_for_3D[0]][2][1]),float(var_3D_binnings[var_for_3D[0]][2][2])],'z':[int(var_3D_binnings[var_for_3D[1]][2][0]),float(var_3D_binnings[var_for_3D[1]][2][1]),float(var_3D_binnings[var_for_3D[1]][2][2])]}
            
            var_for_draw_1 = makeVar(var_3D_formulas[var_for_3D[0]],jetIdxTmp) #jet_pt 
            var_for_draw_2 = makeVar(var_3D_formulas[var_for_3D[1]],jetIdxTmp) #jet_eta

            print 'Plot name: ', plotName
            print 'Axis:      ', axis
            print 'Var x:     ', var
            print 'Var y:     ', var_for_draw_1
            print 'Var z:     ', var_for_draw_2
            
            h_3D = utils.util_funcs.getHisto(v, plotName, axis, [var,var_for_draw_1,var_for_draw_2], cutTmp, sf_final)
            h_3D.Write()

            if r == 'zHFjet' and cat_cut != '':
              for para_var_name,para_var_binning in varBinnings.items():
                for iBin in range(para_var_binning[0]): #0 = nbin
                  binName = para_var_name + '_' + para_var_binning[1][iBin] + '_' + para_var_binning[1][iBin+1]
                  print 'Bin name: ', binName
                  
                  hDir1.cd(binName)
                  plotName = plot + '_' + var_for_3D[0] + '_' + var_for_3D[1] + '_' + binName + '_3D_' + k.replace('_' + chan,'') + '_' + cat_cut.split(':')[0] + '_' + rName
                  para_var_formula = makeVar(varFormulas[para_var_name],jetIdxTmp)
                  cut_bin = '((' + para_var_formula + ' >= ' +  para_var_binning[1][iBin] + ')' + '&&(' + para_var_formula + ' < ' +  para_var_binning[1][iBin+1] + '))'
                  cutTmp_1 = '((' + cutTmp + ')&&' + cut_bin + ')' 
                  print 'Cut bin: ', cut_bin
                  print 'Cut tmp: ', cutTmp_1
                  h_3D = utils.util_funcs.getHisto(v, plotName, axis, [var,var_for_draw_1,var_for_draw_2], cutTmp_1, sf_final)
                  h_3D.Write()

                  hDir1.cd()

          '''      
          #make 1D plots for para_var, for example jet_pt, V_pt
          for para_var_name,para_var_binning in varBinnings.items():
            #In validation region (for example, emu region) only make jet_vtxMass for different jetpt bins
            if r != 'zjet' and r != 'zHFjet' and var != jet_vtxMass_name and para_var_name != var_for_2D[0]: continue #var_for_2D[0] = jet_pt 
            for iBin in range(len(para_var_binning[1])-1):
              plotDir = para_var_name + '_' + para_var_binning[1][iBin] + '_' + para_var_binning[1][iBin+1]
              print 'Plot dir: ', plotDir
              hDir1.cd(plotDir)
              cutTmp_bin = '(' + cutTmp + ')&&(' + + ')'
              
              
              
            
              hDir1.cd()
          '''
    hDir1.cd()
  fOut.cd()

