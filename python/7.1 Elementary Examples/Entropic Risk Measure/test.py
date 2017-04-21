import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import grid
from rpy2.robjects.lib import ggplot2
import numpy as np
import pandas as pd
import rpy2.robjects as robj
import rpy2.robjects.pandas2ri # for dataframe conversion
from rpy2.robjects.packages import importr

# First, make some random data
x = np.random.normal(loc = 5, scale = 2, size = 10)
y = x + np.random.normal(loc = 0, scale = 2, size = 10)
testData = pd.DataFrame( {'x':x, 'y':y} )
print testData

plotFunc = robj.r("""
 library(ggplot2)

function(df){
 p <- ggplot(df, aes(x, y)) +
 geom_point( )

print(p)
 }
""")
gr = importr('grDevices')

# convert the testData to an R dataframe
robj.pandas2ri.activate()
testData_R = robj.conversion.py2ri(testData)

# run the plot function on the dataframe
plotFunc(testData_R)

# ask for input. This requires you to press enter, otherwise the plot
# window closes immediately
#raw_input()
gr.dev_off()

