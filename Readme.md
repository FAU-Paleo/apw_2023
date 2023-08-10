# README

Currently built at [https://fau-paleo.github.io/apw_2022/](https://fau-paleo.github.io/apw_2022/). 

## About Jekyll

The page is built with jekyll, a neet ruby program which translates markdown documents to html and builds these into static websites. It is insanely powerful. The building process itself is exectued by GitHub itself, you make the changes, commit and push them to the repository and in 1 minute they should be visible - given that everything you have done is ok. If you want to learn more about this, check out [this page](https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll). 


If you have Jekyll installed on your computer, you can also render and check out how the page looks like offline - but the page is so simple that this is not necessary.  

## For Instructors 

*You are expected to edit the module-specific pages*, which are in the `_posts` directory. These pages have a front matter (before the dashes), don't worry about those - except for the title, which you are more than welcome to change.

| file                           | correspondent(s)      |
|--------------------------------|-----------------------|
| 2022-08-21-toolset.md          | Adam                  |
| 2022-08-23-paleodiversity.md   | Emma and Wolfgang     |
| 2022-08-25-global-diversity.md | Emma and Adam         |
| 2022-08-27-phylogenetics.md    | Rachel and Laura      |
| 2022-08-29-pyrate.md           | Daniele               |
| 2022-08-31-div-cmr.md          | Lee Hsiang and Isaiah |
| 2022-09-04-paleogeography.md   | Liz and Adam          |
| 2022-09-05-niches.md           | Erin and Tom          |
| 2022-09-07-morphometrics.md    | Ryan                  |
	

- For novice git users: to avoid chaos please only edit the page that you are responsible for and use the GitHub interface. If all of this feels overwhelming, don't worry. I (Adam) is more than happy to put material on for you if you'd like!

- I have included examples on every one of these to 1. include links to external files 2. links to downloadable files and 3. images. You can copy these and modify them to include all the material you want to share. You are also welcome to explore! 

- The data files should be put in the `data` directory, then the format `{{site.baseurl}}/data/<your directory path>/filename.ext` is used to create the link to the file. They need to be committed to the git repo. GitHub has a file size limit, if you want to share files that are bigger than 100MB, please use an external hosting service and provide the download links. I (Adam) can also put up big files here. Do not share sensitive files here. 

- If you have sensitive files that you want to share with the students, please 

- The dates of the posts control when they are visible, future posts are not visible, hence the past dates

- I will make the links to the material in the schedule table available on the fly. To help you explore, they are all live at the moment (filenames are set to a past date last year). Even though the links will not work (not to disorient the students), the URLs will be live. It is a good idea therefore to save the URL of the page that you are editing. 


## Suggested organisation scheme for files

To help keep everything organsied, you'll find subfolders within the main data and slideshows folders under the name of your module(s). We suggest dividing your files as follows:

1. Datasets, scripts, associated files 

These should go in the `data` directory: `{{site.baseurl}}/data/<your directory path>/filename.ext`
The markdown `{{site.baseurl}}/data/1_toolset/metadata.txt` is built as `https://fau-paleo.github.io/apw_2022/data/1_toolset/metadata.txt`

2. Slideshow presentations   
These should go in the `slides` directory: `{{site.baseurl}}/slides/<your directory path>/filename.ext`

The markdown `{{site.baseurl}}/slides/1_toolset/slide_meta.txt` is built as `https://fau-paleo.github.io/apw_2022/slides/1_toolset/slide_meta.txt`


## Markdown

You can find a cheat sheet [here](https://www.markdownguide.org/cheat-sheet/).
