#!/usr/bin/env python
from sys import exit
from runinfo import *
from visit_utils.encoding import encode

# User-defined parameters
db        = "../dumps/dump_*.xmf database"  # relative location of the dump files
var       = "eos3"  # variable name as it appears in the .h5 file
vartext   = "S"  # how you want it to appear in the output
varmin    = 0.0  # min value for color table
varmax    = 30.0  # max value for color table
dumps_pad = 20  # start movie this many steps before bounce
recompute = True  # recompute the individual frames; otherwise, just stitch movie
xmin      =-40.0  # xmin of plot window
xmax      = 40.0  # xmax of plot window
ymin      =-40.0  # ymin of plot window
ymax      = 40.0  # ymax of plot window

# Define expressions
DefineScalarExpression("netheat", "(<Erad0/heat> + <Erad1/heat> + <Erad2/heat>) - (<Erad0/cool> + <Erad1/cool> + <Erad2/cool>)")
DefineScalarExpression("Ftot", "sqrt(<Frad0/dir0/total>^2 + <Frad0/dir1/total>^2)")
DefineVectorExpression("Coords", "coords(mesh)")
DefineScalarExpression("X", "Coords[0]*1.0e5")
DefineScalarExpression("Y", "Coords[1]*1.0e5")
DefineScalarExpression("r", "sqrt(X*X+Y*Y)")
DefineScalarExpression("Ltot", "Ftot*12.566370614359172*r*r")
DefineScalarExpression("vmag", "sqrt(u1*u1+u2*u2)")
DefineVectorExpression("v", "{velocity0,velocity1}")

# DON'T TOUCH BELOW THIS LINE UNLESS YOU KNOW WHAT YOU'RE DOING!!

###############################################################################
# Write color tables
SetCloneWindowOnFirstRef(0)
width, height = 1280, 1024
win = GetGlobalAttributes().windows[GetGlobalAttributes().activeWindow]
ResizeWindow(win, width, height)
SetActiveWindow(win) # Synchronize
size = GetWindowInformation().windowSize
if width < size[0] or height < size[1]:
    ResizeWindow(win, width + (size[0] - width), height + (size[1] - height))
DeleteAllPlots()

# Set save attributes
s = SaveWindowAttributes()
s.format = s.PNG
s.width, s.height = width, height
s.fileName = "rad_ccsn_" + save_text + "_" + vartext + "_"
s.outputToCurrentDirectory = 0
s.outputDirectory = "./frames/"
SetSaveWindowAttributes(s)

# Create plots
# Create plot 1
OpenDatabase(db)
AddPlot("Pseudocolor", var, 0, 0)
atts = PseudocolorAttributes()
atts.colorTableName = "Spectral"
atts.invertColorTable = 1
atts.centering = atts.Nodal
atts.minFlag = 1
atts.min     = varmin
atts.maxFlag = 1
atts.max     = varmax
SetPlotOptions(atts)
AddOperator("Reflect", 0)
opAtts = ReflectAttributes()
opAtts.reflections = (1,1,0,0,0,0,0,0)
SetOperatorOptions(opAtts)
silr = SILRestriction()
silr.TurnOnAll()
SetPlotSILRestriction(silr, 0)

SetActivePlots(0)

DrawPlots()

# Set the view
view = View2DAttributes()
view.windowCoords = (-250, 250, -250, 250)
view.viewportCoords = (0.2, 0.95, 0.1, 0.9)
view.fullFrameActivationMode = view.Auto  # On, Off, Auto
view.fullFrameAutoThreshold = 100
view.xScale = view.LINEAR  # LINEAR, LOG
view.yScale = view.LINEAR  # LINEAR, LOG
view.windowValid = 1
SetView2D(view)

# Set the general annotation attributes
annot = AnnotationAttributes()
annot.userInfoFlag = 0
annot.databaseInfoTimeOffset = -tbounce
annot.axes2D.autoSetScaling = 0
annot.axes2D.xAxis.title.userUnits = 1
annot.axes2D.xAxis.title.units = "km"
annot.axes2D.yAxis.title.userUnits = 1
annot.axes2D.yAxis.title.units = "km"
SetAnnotationAttributes(annot)

# Set text annotation (title) attributes
win0_obj002 = CreateAnnotationObject("Text2D", "Title")
win0_obj002.visible = 1
win0_obj002.position = (0.35, 0.91)
win0_obj002.height = 0.03
win0_obj002.textColor = (0, 0, 0, 255)
win0_obj002.useForegroundForTextColor = 1
win0_obj002.text = title_text
win0_obj002.fontFamily = win0_obj002.Arial  # Arial, Courier, Times
win0_obj002.fontBold = 0
win0_obj002.fontItalic = 0
win0_obj002.fontShadow = 0

# Set legend (colorbar) attributes
win0_legend000 = GetAnnotationObject(GetPlotList().GetPlots(0).plotName)
win0_legend000.active = 0
win0_legend000.managePosition = 1
win0_legend000.position = (0.05, 0.9)
win0_legend000.xScale = 1
win0_legend000.yScale = 1
win0_legend000.textColor = (0, 0, 0, 255)
win0_legend000.useForegroundForTextColor = 1
win0_legend000.drawBoundingBox = 0
win0_legend000.boundingBoxColor = (0, 0, 0, 50)
win0_legend000.numberFormat = "%# -9.4g"
win0_legend000.fontFamily = win0_legend000.Arial  # Arial, Courier, Times
win0_legend000.fontBold = 0
win0_legend000.fontItalic = 0
win0_legend000.fontShadow = 0
win0_legend000.fontHeight = 0.015
win0_legend000.drawLabels = win0_legend000.Values # None, Values, Labels, Both
win0_legend000.drawTitle = 1
win0_legend000.drawMinMax = 1
win0_legend000.orientation = win0_legend000.VerticalRight  # VerticalRight, VerticalLeft, HorizontalTop, HorizontalBottom
win0_legend000.controlTicks = 1
win0_legend000.numTicks = 5
win0_legend000.minMaxInclusive = 1
win0_legend000.suppliedValues = ()
win0_legend000.suppliedLabels = ()

SetActiveWindow(GetGlobalAttributes().windows[0])

if recompute:
    # Cycle through states and save each window
    for state in range(dump_start-dumps_pad,dump_stop+1):
        SetTimeSliderState(state)
        SaveWindow()

encode('frames/'+s.fileName+'%04d.png',s.fileName+'movie.mp4')

ClearCacheForAllEngines()
DeleteAllPlots()
CloseDatabase(db)
ClearCacheForAllEngines()
CloseComputeEngine()
sys.exit()

