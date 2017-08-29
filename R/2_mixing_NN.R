# setup ----------------

library(tensorflow)

#flags <- tf$app$flags
#flags$DEFINE_integer('n_hidden_layer', 5L, 'Number of hidden layers.')
#flags$DEFINE_integer('n_nodes_per_layer', 10L, 'Number of nodes per hidden layer.')
#f#lags$DEFINE_float('p_keep', 0.8, 'Expecatation of fration of not droppedout units.')
#flags$DEFINE_float('f_vali', 1/5, 'fraction of data used for validation')
#flags$DEFINE_integer('batch_size', 100, 'batch size')
#FLAGS <- parse_flags()
FLAGS <- list()
FLAGS$n_hidden_layer<-1L;
FLAGS$n_nodes_per_layer<-100L;
FLAGS$p_keep<-0.8;
FLAGS$f_vali<-1/5;
FLAGS$batch_size<-1000;

tf$logging$set_verbosity(tf$logging$INFO)


tRunif <- tf$random_uniform(shape(length(E),1))
batch <-
tX <- tf$constant(E, tf$float32)# tf$placeholder(tf$float32, shape(NULL))
tTargets <- tf$constant(as.numeric(y), tf$float32)# tf$placeholder(tf$float32, shape(NULL))


net <- tRunif

for(i in 1:FLAGS$n_hidden_layer) {
  net = tf$contrib$layers$fully_connected(
    net, FLAGS$n_nodes_per_layer, activation_fn=tf$nn$elu,
    weights_initializer=tf$contrib$layers$xavier_initializer())
  net = tf$nn$dropout(net, keep_prob=FLAGS$p_keep)
}
net = tf$contrib$layers$fully_connected(net, 1L, activation_fn=NULL,weights_initializer=tf$contrib$layers$xavier_initializer())/10

loss=tf$tanh(tf$reduce_mean(tf$nn$log_poisson_loss(targets = tTargets, log_input =  tf$squeeze(net,axis=1L)-tf$log(tX),  compute_full_loss=FALSE)))

theta = tf$exp(tf$squeeze(net,axis=1L))
d_theta_d_runif <- tf$squeeze(tf$gradients(theta, tRunif))
density = tf$div(1, d_theta_d_runif)

lambda <- 1
monotonicty_loss = lambda*tf$tanh(tf$reduce_mean(1/(1+tf$abs(d_theta_d_runif))) + tf$reduce_mean(tf$abs(d_theta_d_runif)- d_theta_d_runif))

optimizer <- tf$train$AdamOptimizer(learning_rate=0.01)
train <- optimizer$minimize(loss+monotonicty_loss)

# Launch the graph and initialize the variables.
config = tf$ConfigProto(
  device_count = dict('GPU'=0L)
)
sess = tf$Session(config=config)
#sess = tf$Session()
sess$run(tf$global_variables_initializer())
i <- 0
loss_log <- sess$run(loss)
monotonicty_loss_log <- sess$run(monotonicty_loss)
# training -----------------
for ( j in seq_len(200)) {
  i <- i+1
  res <- sess$run(list(train, loss, monotonicty_loss))
  loss_log[i+1] <- res[[2]]
  monotonicty_loss_log[i+1] <- res[[3]]
  if (i %% 50 == 0)
    plot(rep(seq(0, i),2),c(loss_log[1:(i+1)],monotonicty_loss_log[1:(i+1)]), pch=21, bg=rep(c("red","green3"),each=i+1))
}
# -------------------------------

sess$run(tf$exp(tf$squeeze(net,axis=1L)-tf$log(tX)))
sess$run(loss)
sess$run(loss)
res <- sess$run(list(theta,density, tRunif))
plot(res[[1]], res[[2]], log = "y")
plot(res[[3]], res[[1]])
saver <- tf$train$Saver()
saver$save(sess,'test_w/g')

test_writer = tf$summary$FileWriter('test_w')
test_writer$add_graph(net$graph)


# -------------------------------
sess$close()
