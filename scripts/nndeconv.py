from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import tensorflow as tf
import time
from bunch import Bunch
from sklearn.model_selection import train_test_split

tf.logging.set_verbosity(tf.logging.INFO)


FLAGS = Bunch()
FLAGS.batch_size = 10
FLAGS.p_use = 1.0 #.9
FLAGS.n_nodes_per_layer = 1000
FLAGS.n_epochs = 100
FLAGS.test_every = 100
FLAGS.train_dir = "nndeconv_out"
FLAGS.data_file = "../data_raw/deconv.csv"
FLAGS.lambda_monotonicity = 0.1
FLAGS.learning_rate = 0.02
FLAGS.lambda_L2 = 0.01
FLAGS.n_int = 100

# Load datasets.
data = np.loadtxt(FLAGS.data_file, skiprows=1, delimiter=",", dtype=np.float32)
train, test = train_test_split(data, test_size=0.2)
E_train = tf.constant(train[:, 1], name="exposure")
y_train = tf.constant(train[:, [2]], name="counts")
E_test = tf.constant(test[:, 1], name="exposure")
y_test = tf.constant(test[:, [2]], name="counts")

[E_train_batch, y_train_batch] = tf.train.batch(tf.train.slice_input_producer([E_train, y_train], num_epochs=FLAGS.n_epochs),
                                                batch_size=FLAGS.batch_size)
# [E_test_batch, y_test_batch] = tf.train.batch(tf.train.slice_input_producer([E_test, y_test]), batch_size=FLAGS.batch_size)


modes = ["train","test","predict","sample"]

L2_loss = 0

W1 = tf.get_variable("W1", shape=[1, FLAGS.n_nodes_per_layer ], initializer=tf.contrib.layers.xavier_initializer())
tf.summary.histogram("W1", W1, collections=modes)
L2_loss += tf.norm(W1)
b1 = tf.get_variable("b1", shape=[FLAGS.n_nodes_per_layer], initializer=tf.zeros_initializer())
tf.summary.histogram("b1", b1, collections=modes)
W2 = tf.get_variable("W2", shape=[FLAGS.n_nodes_per_layer , 1], initializer=tf.contrib.layers.xavier_initializer())
L2_loss += tf.norm(W2)
L2_loss *= FLAGS.lambda_L2


tf.summary.scalar('l2_loss', L2_loss, collections=modes)
b2 = tf.get_variable("b2", shape=[1], initializer=tf.zeros_initializer())
tf.summary.histogram("W2", W2, collections=modes)
tf.summary.histogram("b2", b2, collections=modes)

E = {"train": E_train_batch,
     "test": E_test,
     "predict": tf.placeholder(shape=[None], dtype=tf.float32)}

y = {"train": y_train_batch,
     "test": y_test,
     "predict": tf.placeholder(shape=[None, 1], dtype=tf.float32)}
p_use = {"train": FLAGS.p_use,
         "test": 1,
         "predict": tf.placeholder_with_default(np.float32(1), shape=[]),
         "sample": 1}

n_sample = tf.placeholder(dtype=tf.int32, shape=[None], name="n_sample")

loss_data = dict()
monotonicity_loss = dict()
summaries = dict()
for m in modes:
    R = tf.random_uniform(shape=tf.concat([[FLAGS.n_int], tf.shape(y[m])], axis=0) if m is not "sample" else n_sample)
    net = R
    tf.summary.histogram("in"+m, net, collections=[m])
    net = tf.nn.elu(tf.add(tf.tensordot(net, W1, axes=1), b1))
    tf.summary.histogram("hidden_layer"+m, net, collections=[m])
    if p_use[m] is not 1:
        net = tf.nn.dropout(net, keep_prob=p_use[m])
        tf.summary.histogram("hidden_layer_after_dropout"+m, net, collections=[m])
    net = tf.add(tf.tensordot(net, W2, axes=1), b2)

    theta = tf.exp(net)
    d_theta_d_runif = tf.gradients(theta, R)
    density = tf.divide(1, tf.abs(d_theta_d_runif))
    if m == "sample":
        sample_R = tf.identity(R, "sample_R")
        sample_theta = tf.identity(theta, "sample_theta")
        sample_density = tf.identity(density, "sample_density")
    else:
        tf.summary.histogram("out" + m, net, collections=[m])
        loss_data[m] = tf.reduce_mean(-tf.log(tf.reduce_mean(tf.exp(-tf.nn.log_poisson_loss(
            targets=y[m],
            log_input=net + tf.expand_dims(tf.expand_dims(tf.log(E[m]), 1), 0),
            compute_full_loss=True)), axis=0)))
        tf.summary.scalar('loss_data' + m, loss_data[m], collections=[m])
        monotonicity_loss[m] = FLAGS.lambda_monotonicity * tf.reduce_mean(tf.abs(d_theta_d_runif) - d_theta_d_runif)
        tf.summary.scalar('monotonicity_loss'+m, monotonicity_loss[m], collections=[m])
        tf.summary.scalar('total_loss'+m, loss_data[m] + monotonicity_loss[m]+ L2_loss, collections=[m])
    summaries[m] = tf.summary.merge_all(key=m)

optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)
train = optimizer.minimize(loss_data["train"] + monotonicity_loss["train"]+ L2_loss)

# Create a saver for writing training checkpoints.
saver = tf.train.Saver()

# Create the op for initializing variables.
init_op = tf.group(tf.global_variables_initializer(),
                   tf.local_variables_initializer())

# Launch the graph and initialize the variables.
config = tf.ConfigProto(device_count={'GPU': 0})
with tf.Session(config=config)  as sess:
    sess.run(init_op)
    writer = tf.summary.FileWriter(FLAGS.train_dir + "/"+time.strftime("%Y-%m-%d %H-%M-%S", time.gmtime()))

    # Start input enqueue threads.
    coord = tf.train.Coordinator()
    threads = tf.train.start_queue_runners(sess=sess, coord=coord)
    i = 0
    writer.add_summary(sess.run(summaries["train"]), i)
    writer.add_graph(tf.get_default_graph())
    writer.add_summary(sess.run(summaries["test"]), i)
    # training -----------------
    try:
        while not coord.should_stop():
            i = i + 1
            [_, sum] = sess.run([train, summaries["train"]])
            writer.add_summary(sum, i)
            if i % FLAGS.test_every is 0:
                writer.add_summary(sess.run(summaries["test"]), i)
                print(i)
            # -------------------------------
            # Save a checkpoint periodically.
            if i % 1000 == 0:
                print('Saving')
                saver.save(sess, FLAGS.train_dir+"/out.sav", global_step=i)
    except tf.errors.OutOfRangeError:
        print('Saving')
        saver.save(sess, FLAGS.train_dir+"/final.sav", global_step=i)
        print('Done training for %d epochs, %d steps.' % (FLAGS.n_epochs, i))
    finally:
      # When done, ask the threads to stop.
      coord.request_stop()

    # Wait for threads to finish.
    coord.join(threads)