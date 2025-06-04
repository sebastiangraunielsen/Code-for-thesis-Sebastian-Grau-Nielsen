import numpy.random as nr
import random
import numpy.linalg as alg
import scipy.linalg as salg
import tensorflow as tf
from tensorflow.keras import layers, losses, models
from tensorflow.keras.models import Model
from tensorflow.keras.initializers import TruncatedNormal, Constant
from tensorflow.keras.callbacks import Callback
import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import cmath
from PIL import Image
import idx2numpy
import os

n=20 # number of training samples
phalf=28 # number of pixels in each column or each row
p=phalf*phalf # number of pixels (phalf x phalf)
epochs=2000 # maximum number of epochs for learning CNN
std_noise = 0.5 # standard deviation for Gaussian noise

# Read the data for 3x3 CNN
# Load MNIST dataset
(x_train, y_train), (_, _) = tf.keras.datasets.mnist.load_data()

# Normalize to [0,1]
x_train = x_train.astype('float32') / 255.

# Add Gaussian noise
x_train_noisy = x_train + std_noise * np.random.normal(loc=0.0, scale=1.0, size=x_train.shape)

# Clip values to [0,1]
x_train_noisy = np.clip(x_train_noisy, 0., 1.)

# Now, find the images to train on
index=[[] for i in range(10)]

for i in range(10):
  index[i]=list(np.where(y_train==i))[0]

x_train_reduced = np.zeros((n, phalf, phalf))
x_train_noisy_reduced = np.zeros((n, 28, 28))
x_test_reduced = np.zeros((100, 28, 28))
x_test_noisy_reduced = np.zeros((100, 28, 28))
  
for i in range(n):
  x_train_reduced[i] = x_train[index[i%10][i],:,:].copy()
  x_train_noisy_reduced[i] = x_train_noisy[index[i%10][i],:,:].copy()

for i in range(100):
  x_test_reduced[i] = x_train[index[i%10][i+n],:,:].copy()
  x_test_noisy_reduced[i] = x_train_noisy[index[i%10][i+n],:,:].copy()

x_train_reduced = np.reshape(x_train_reduced, (n, 28, 28, 1))
x_train_noisy_reduced = np.reshape(x_train_noisy_reduced, (n, 28, 28, 1))
x_test_reduced = np.reshape(x_test_reduced, (100, 28, 28, 1))
x_test_noisy_reduced = np.reshape(x_test_noisy_reduced, (100, 28, 28, 1))

model_3layer = models.Sequential([
    layers.Input(shape=(28, 28, 1)),
    
    layers.Conv2D(32, (3, 3), padding='same'),
    layers.Activation('relu'),

    layers.Conv2D(64, (3, 3), padding='same'),
    layers.Activation('relu'),

    layers.Conv2D(1, (3, 3), padding='same'),
    layers.Activation('sigmoid')
])


model_28 = models.Sequential([
  layers.Input(shape=(28, 28, 1)),
  layers.Conv2D(32, (28, 28), padding='same'),
  layers.Activation('relu'),
  
  layers.Conv2D(64, (28, 28), activation='relu', padding='same'),
  layers.Activation('relu'),
  
  layers.Conv2D(1, (3,3), activation='sigmoid', padding='same'),
  layers.Activation('sigmoid')

 ])



# Compile the models
model_3layer.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.0001), loss='mean_squared_error')
model_28.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.0001), loss='mean_squared_error')

# Train the model
errlist_3layer = []
errlist_28 = []

for epoch in range(1, epochs + 1):
    # Train for one epoch
    model_3layer.fit(x_train_noisy_reduced, x_train_reduced, epochs=1, batch_size=64, verbose=0)
    model_28.fit(x_train_noisy_reduced, x_train_reduced, epochs=1, batch_size=64, verbose=0)

    # Every 10th epoch, compute mean Frobenius norm on test set
    if epoch % 10 == 0:
        y_pred_3 = model_3layer.predict(x_test_noisy_reduced, verbose=0)
        y_pred_28 = model_28.predict(x_test_noisy_reduced, verbose=0)
        frob_norm_3layer = 0
        frob_norm_28 = 0

        # Compute Frobenius norm for each test image
        for i in range(100):
          frob_norm_3layer = frob_norm_3layer + alg.norm(y_pred_3[i].reshape([p]) - x_test_reduced[i].reshape([p]))
          frob_norm_28 = frob_norm_28 + alg.norm(y_pred_28[i].reshape([p])-x_test_reduced[i].reshape([p]))

        errlist_3layer = np.append(errlist_3layer, np.array([frob_norm_3layer])/100)
        errlist_28 = np.append(errlist_28, np.array([frob_norm_28])/100)

epochs_range = np.reshape(np.array([range(10, epochs + 1, 10)]), (len(errlist_3layer)))
plt.plot(epochs_range, errlist_3layer, label='3-layer 3x3 CNN')
plt.plot(epochs_range, errlist_28, label='2-layer 28x28 + 1-layer 3x3 CNN')
plt.xlabel('Epochs')
plt.ylabel('Mean test error')
plt.show
