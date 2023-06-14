import os

PATH = 'figure'

def get_img(path):
    img_list = []
    for dirs, subdirs, files in os.walk(path):
        for f in files:
            img_list.append(os.path.join(dirs, f))
    return img_list

# remove all files in folders
img_lst = get_img(PATH)
for img in img_lst:
    os.remove(img)