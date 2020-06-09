# Imageanalysis
# Final Project on Counting Puncta in Mouse Colonocytes
## For my final project, I will be observing colocalization of puncta in Rab proteins with and without alpha toxins. The Red, Green and Blue chnnaels have different information: 
           R: Rab protein marker
           G: Alpha tocxin marker
           B: Nucleus marker
## Introduction 
  ### Investigating a code that will allow a scientist to gather puncta from mouse colonocytes with the treated mice with alpha toxin and without alpha toxin. These mice will also have the Rab proteins 5 and 11. The goal is to decipher whether or not there is colocalization in Rab proteins with alpha toxins and Rab proteins within the nucleus. 

## How to start
### 1) Import all packages necessary for extracting files and counting puncta
      %matplotlib inline
      from matplotlib import pyplot as plt
      import numpy as np

      from PIL import Image
      from skimage.feature import blob_dog, blob_log, blob_doh
      from skimage.color import rgb2gray
      from skimage.draw import circle

      import os
      import glob

### 2) Create code to investigate the colocalization of the Rab protein with and without alpha toxin puncta. Calling out the images from your folder
      import imageio
      import glob

      e_im = []
      c_im = []
      for image_path in     glob.glob("/Users/juantirado/Desktop/Progamming/Final_Project_Folder/Rab5/Experimental/Rab5_alpha_toxin3_not_counted.png"):
          e_im = imageio.imread(image_path)

      for image_path in glob.glob("/Users/juantirado/Desktop/Progamming/Final_Project_Folder/Rab5/Control/Rab5_no_toxin3_not_counted.png"):
          c_im = imageio.imread(image_path)
          
### 3) Importing the image using this call function for the control image
      plt.imshow(c_im)
      print("x,y,RGB - >", c_im.shape)
### 4) Change image to grayscale from the control image
      image_gray = rgb2gray(c_im)
      plt.imshow(image_gray,cmap="gray")
### 5) Gathering data from the control file
    def plot_blobs(img,blobs):
    """
    Plot a set of blobs on an image.
    """
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 

    ax.imshow(img, interpolation='nearest')
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color="yellow", linewidth=2, fill=False)
        ax.add_patch(c)
### 6) Utilizing the blob detection techniques with dog,log,and doh for the control image. This function also allows for a display of 3 images detecting the blobs
    blobs_log = blob_log(rgb2gray(c_im), max_sigma=30, num_sigma=10, threshold=.1)
    blobs_log[:, 2] = blobs_log[:, 2] * np.sqrt(2)

    blobs_dog = blob_dog(rgb2gray(c_im), min_sigma = 20, max_sigma = 50, threshold=.3)
    blobs_dog[:, 2] = blobs_dog[:, 2] * np.sqrt(2)

    blobs_doh = blob_doh(rgb2gray(c_im), max_sigma=30, threshold=.01)
    blobs_list = [blobs_log, blobs_dog, blobs_doh]
    colors = ['yellow', 'white', 'blue']
    titles = ['Laplacian of Guassian', 'Difference of Gaussian', 
              'Determinant of Hessian']

    sequence = zip(blobs_list, colors, titles)

    fig, axes =plt.subplots(1,3,figsize=(9,3), sharex=True, sharey=True)
    ax = axes.ravel()

    for idx, (blobs, color, title) in enumerate(sequence):
        ax[idx].set_title(title)
        ax[idx].imshow(c_im)
        for blob in blobs:
            y, x, r =blob
            c = plt.Circle((x,y), r, color=color, linewidth=2, fill=False)
            ax[idx].add_patch(c)
        ax[idx].set_axis_off()

    plt.tight_layout()
    plt.show()
### 7) This code function gathers 5 nuclei from blobs dog from the control group and the following confirms the gathering of 5 nuclei using a tuple.
    nuclei = []
    for counter, item in enumerate(blobs_dog):
        contender = blobs_dog[np.random.randint(2, len(blobs_dog))]
        if len(nuclei) < 5 and contender[1] >35:
            nuclei.append(contender)
            
     nuclei = tuple(nuclei)
     print(len(nuclei))
### 8) This code allows for detection of 5 nuclei 
      blobs_list = [nuclei]
      colors = ['white']
      titles = ['Difference of Gaussian']
      sequence = zip(blobs_list, colors, titles)
      blobs_list=nuclei


      fig, axes =plt.subplots(1,3,figsize=(9,3), sharex=True, sharey=True)
      ax = axes.ravel()

      for idx, (blobs, color, title) in enumerate(sequence):
          ax[idx].set_title(title)
          ax[idx].imshow(c_im)
          for blob in blobs:
              y, x, r =blob
              c = plt.Circle((x,y), r, color=color, linewidth=2, fill=False)
              ax[idx].add_patch(c)
          ax[idx].set_axis_off()

      plt.tight_layout()
      plt.show()
### 9) This code allows a preparation of a circle mask for each sampled nuclei and storing in a dictionary
      circle_masks = {}
      for counter, nucleus in enumerate(nuclei):
          circle_masks.setdefault(counter, circle(nucleus[0], nucleus[1], nucleus[2]))
### 10) This code allows for gathering a blob color function 
      def get_blob_color(image,blob):

      print(blob)

      # Grab circle center (a,b) and radius (r)
      (a, b, r) = blob

      # Draw circle for blob
      circle_mask = circle(a,b,r)

      # Create mask of False over whole image
      mask = np.zeros(image.shape[0:2],dtype=np.bool)
      mask[circle_mask] = True

      num_pixels = np.sum(mask)
      red = np.sum(image[mask,0])/num_pixels
      green = np.sum(image[mask,1])/num_pixels
      blue = np.sum(image[mask,2])/num_pixels

      return red, green, blue
 ### 11) This code allows to obtain specific colored puncta. In this case we are gathering red puncta because there is no alpha toxin (green channel)
    red_puncta = []
    not_red = 0
    red = 0
    for counter, puncta in enumerate(nuclei):
        if blobs[counter][0] < 100:
            not_red += 1
        else:
            red += 1
            red_puncta.append(puncta)

    print(f'There are {red} puncta that are red and {not_red} false positives.')
### 12) Storing the puncta in a dictionary
      orted_control_puncta = {}
      for counter, nucleus in enumerate(nuclei):
          center_x = nucleus[0]
          center_y = nucleus[1]
          radius = nucleus[2]
          contained_puncta = []
          for _, puncta in enumerate(nuclei):
              x = puncta[0]
              y = puncta[1]      
              if (x-center_x)**2 + (y-center_y)**2 <= radius**2:
                  contained_puncta.append(puncta)
          sorted_control_puncta.setdefault(counter, contained_puncta)

      #### Verify that all puncta were sorted
      verification = np.empty(len(nuclei))
      for counter, nucleus in enumerate(nuclei):
          verification[counter] = len(sorted_control_puncta[counter])

      if np.sum(verification) == len(nuclei):
          print(f'All {len(nuclei)} puncta were successfully sorted.')
      else:
          print('At least one of the puncta was not properly sorted.')
### 13) Running data analysis on  the control group of Rab5 protein
      puncta_per_nucleus_control = [len(nucleus) for nucleus in sorted_control_puncta.values()]
      print(puncta_per_nucleus_control)
      
      plt.hist(puncta_per_nucleus_experimental)
      
      experimental_mean = np.mean(puncta_per_nucleus_experimental)
      print(experimental_mean)
      
      
## Repeat the process for any Rab protein and image and confirm your analysis




