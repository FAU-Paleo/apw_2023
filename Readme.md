# README

Currently built at [https://fau-paleo.github.io/apw_2023/](https://fau-paleo.github.io/apw_2023/). 

## About Jekyll

The page is built with jekyll, a very neat ruby program which translates markdown documents to html and builds these into static websites. It is insanely powerful. The building process is executed by GitHub itself, you make the changes, commit and push them to the repository and in 1 minute they should be visible - given that everything you have done is ok. If you want to learn more about this, check out [this page](https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll). 

If you have Jekyll installed on your computer, you can also render and check out how the page looks offline - but the page is so simple that this is not necessary.  

## Sharing material

All actual files of the webpage are in the `docs` directory. **You are expected to edit the module-specific pages**, which are in the `docs/_posts` directory. These pages have a front matter (before the dashes), don't worry about those - except for the title, which you are welcome to change. Here is a table of what files were supposed to be edited by whom:

| file                         | correspondent(s)        | link                                                                                                         |
|------------------------------|-------------------------|--------------------------------------------------------------------------------------------------------------|
| 2022-08-21-toolset.md        | Adam, Emma              | [https://fau-paleo.github.io/apw_2023/toolset/](https://fau-paleo.github.io/apw_2023/toolset/)               |
| 2022-08-23-paleodiversity.md | Emma, Wolfgang and Adam | [https://fau-paleo.github.io/apw_2023/paleodiversity/](https://fau-paleo.github.io/apw_2023/paleodiversity/) |
| 2022-08-27-phylogenetics.md  | Rachel                  | [https://fau-paleo.github.io/apw_2023/phylogenetics/](https://fau-paleo.github.io/apw_2023/phylogenetics/)   |
| 2022-08-29-pyrate.md         | Daniele                 | [https://fau-paleo.github.io/apw_2023/pyrate_ai/](https://fau-paleo.github.io/apw_2023/pyrate_ai/)           |
| 2022-08-31-div-cmr.md        | Lee Hsiang and Isaiah   | [https://fau-paleo.github.io/apw_2023/cmr/](https://fau-paleo.github.io/apw_2023/cmr/)                       |
| 2022-09-04-paleogeography.md | Liz and Adam            | [https://fau-paleo.github.io/apw_2023/paleogeography/](https://fau-paleo.github.io/apw_2023/paleogeography/) |
| 2022-09-05-niches.md         | Erin and Tom            | [https://fau-paleo.github.io/apw_2023/niches/](https://fau-paleo.github.io/apw_2023/niches/)                 |
| 2022-09-07-morphometrics.md  | Ryan                    | [https://fau-paleo.github.io/apw_2023/morphometrics/](https://fau-paleo.github.io/apw_2023/morphometrics/)   |


- Material from the previous year are avaiable at this GitHub repo: [GitHub previous year](https://github.com/FAU-Paleo/apw_2022). You can copy markdown from there if you wish. For the sake of convenience I included the old files (`images`, `data` and `slides` directories) in the new repo which are in the `docs/_legacy` directory.

- For novice git users: to avoid chaos please only edit the page that you are responsible for and use the GitHub interface. If all of this feels overwhelming, don't worry. I (Adam) am happy to put the material on the webpage for you if you need help!

- I have included examples in every one of these files to 1. include links to external files 2. links to downloadable files and 3. images. You can copy these and modify them to include all the material you want to share. You are also welcome to explore the first example `docs/_posts/2022-08-21-toolset.md` where I have started to put things already from the previous year.  

- The data files should be put in the `data` directory, then the format `{{site.baseurl}}/data/<your directory path>/filename.ext` is used to create the link to the file. The `{{site.baseurl}}` bit will become `/apw_2023`. The data files also need to be committed to the git repo. GitHub has a file size limit, if you want to share files that are bigger than 100MB, please use an external hosting service and provide the download links. I (Adam) can also put up big files here. 

- If you have sensitive files that you want to share with the students, please use private, external links to these (e.g. dropbox/google drive). Do not add them to the repo because that is relatively difficult to delete. Once the course is over, links to such sensitive files can be invalidated and nobody else will be able to access them.

- The dates of the posts control when they are visible, future posts are not visible, hence the past dates. I will make the links to the material in the schedule table available on the fly. To help you explore, they are all live at the moment (filenames are set to a past date last year). Even though the links will not work (not to disorient the students), the URLs will be live. It is a good idea therefore to save the URL of the page that you are editing. 


## Suggested organisation scheme for files

To help keep everything organsied, you'll find subfolders within the main data directories under the name of your module(s). These should go in the `docs/data` directory: `{{site.baseurl}}/data/<your directory path>/filename.ext`
The file `{{site.baseurl}}/data/1_toolset/metadata.txt` is built as `https://fau-paleo.github.io/apw_2023/data/1_toolset/metadata.txt`


## Markdown

You can find a cheat sheet [here](https://www.markdownguide.org/cheat-sheet/).
