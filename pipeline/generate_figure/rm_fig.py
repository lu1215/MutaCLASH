import os

PATH = ['figure','log']

def get_img(path):
    img_list = []
    for dirs, subdirs, files in os.walk(path):
        for f in files:
            img_list.append(os.path.join(dirs, f))
    return img_list

# remove all files in folders
for folder in PATH:
    img_lst = get_img(folder)
    for img in img_lst:
        os.remove(img)