# A-Two-Stage-Machine-Learning-Framework-for-Predicting-Sporadic-E-Occurrence-and-Intensity
“A Two-Stage Machine Learning Framework for Predicting Sporadic E Occurrence and Intensity“ code
The framework consists of two sequential neural networks designed to mitigate the "dilution effect" caused by the overwhelming number of non-Es samples.
Stage 1 net_occ: A fully connected network estimating Es occurrence probability P. It includes a Batch Normalization layer and three hidden linear layers (42 units each) with ReLU activations, followed by a Sigmoid output.
Stage 2 net_reg: A conditional regression network activated only for positive events via a Gate Function (P > Threshold). It employs three hidden linear layers (52 units each) to predict conditional intensity Y_cond.
