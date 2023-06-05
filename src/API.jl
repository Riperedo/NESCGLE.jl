include("NonInstantaneousProcess/CoolingRate/API.jl")
include("NonInstantaneousProcess/Hysteresis/API.jl")

logo = "\n     0xlccoxkO0K        
   kl;''''';okkkOK      #  ### #  # # #   # ### ###
 Oc,''''''''l00kkkkx0   #  # # ## # # ## ## #   #
k;''''''',:oK  XOkkl:O  #  ### # ## # # # # ##  ##
:'''',cdk0X     Kkx:'cK ## # # #  # # #   # #   ###
''',ck          Xkl,',k 
'',lOX          kc,'';k Laboratorio Nacional de la
:'cxkK     X0Odc,''''cK Ingeniería de la Materia
k:lkkOX  0o:,''''''':O  Fuera del Equilibrio
 Oxkkkk00c'''''''',l0   
   XKOkkko;''''';lO     HandsOn NE-SCGLE Version 0.02
      KOOkxlccld0       (11-05-2023)"

CodeMap = "NE-SCGLE
  │
  ├─SCGLE
  │  │
  │  ├─StructureFactor
  │  │  │
  │  │  ├─HardSphere
  │  │  │  │
  │  │  │  │-PercusYevick
  │  │  │  │-VerletWeis
  │  │  │  │-BlipFunction
  │  │  -  │-RandomPhaseApproximation
  │  │  │-Mixtures
  │  │  │-DipolarHardSpheres
  │  │
  │  ├─Grids
  │  │
  │  ├─Dynamics
  │  │  │-Steps
  │  │  │-TransportProperties
  │  │
  │  └─ArrestDiagram
  │     |-ArrestLine
  │     -
  │     |-Spinodal
  │
  ├─InstataneousProcess
  │  │-AdaptativeGrid
  │  │-NETP
  │
  └─NonInstatanoeusProcess
     |-CoolingRate
     -
     |-Hysteresis"
