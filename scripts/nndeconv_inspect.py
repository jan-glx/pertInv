
import tensorflow as tf
from bunch import Bunch

FLAGS = Bunch()
FLAGS.train_dir = "nndeconv_out"


sess=tf.Session()
#First let's load meta graph and restore weights
saver = tf.train.import_meta_graph(FLAGS.train_dir+'/final.sav-16000.meta')
saver.restore(sess,tf.train.latest_checkpoint(FLAGS.train_dir))

graph = tf.get_default_graph()

sample_R = graph.get_tensor_by_name("sample_R:0")
sample_theta = graph.get_tensor_by_name("sample_theta:0")
sample_density = graph.get_tensor_by_name("sample_density:0")
n_sample = graph.get_tensor_by_name("n_sample:0")

[r,t,d] = sess.run([sample_R,sample_theta,sample_density], feed_dict={n_sample: [10000,1]})

import matplotlib.pyplot as plt


# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].scatter(t, r)
axarr[0].set_ylim([0, 1])
axarr[0].set_title('Sharing X axis')
axarr[1].scatter(t, d)

