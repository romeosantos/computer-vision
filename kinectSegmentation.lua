require 'Kinect'
require 'xLearn'
mincut = require 'mincut'
k = Kinect(640,480)
k:getRGBD()
image.display(k.rgb)
output = torch.Tensor(640,480,3)
mincut.segmentation(k.rgb,k.depth,output)
image.display(output)