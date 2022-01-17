# microglia-phenotyping

SOP Micro glia activation:

1)	Create your own directory (e.g. your name) and copy in the recent Omnisphero folder (Box…Lein Lab shared files/Software/Omnisphero

2)	Open Matlab2014b

3)	Close all open script windows

4)	Go to open (or Ctrl+O) 

![Picture1](https://user-images.githubusercontent.com/97848410/149844249-85e12489-3779-450d-aa42-84c332d54cdc.jpg)

5)	Browse to your Omnisphero directory and select the Omnisphero file of Type MATLAB Code

![Picture2](https://user-images.githubusercontent.com/97848410/149844259-b405a5bb-7630-4d2f-a396-dcda0aeff27d.jpg)

6)	A green script window will appear

![Picture3](https://user-images.githubusercontent.com/97848410/149844261-693f9e25-54b6-4aef-9d23-4ea9c8e4fd63.jpg)

7)	Hit run (green triangle)

8)	Press Change folder

![Picture4](https://user-images.githubusercontent.com/97848410/149844263-2970e2a9-835c-4ac3-bfc8-aa55bbe4a2af.jpg)

9)	The GUI of Omnisphero will pop up

10)	Preprocess the images now by going to LoadMicrogliaMorphology

![Picture5](https://user-images.githubusercontent.com/97848410/149844264-7dee26c0-a0f3-4eff-9e7a-8b43f5ef04b8.jpg)

a.	The following question will pop up. Chose no for fluorescent images

![Picture6](https://user-images.githubusercontent.com/97848410/149844266-b02f895c-e121-447b-b7c5-5828d272f048.jpg)

b.	Fluorescent images Two questions will pop up

![Picture7](https://user-images.githubusercontent.com/97848410/149844267-629bcb2b-acc8-4e84-819f-d6870be86293.jpg)

Hit yes if you have CD68 fluorescence image (yes in case of sample data provided)

![Picture8](https://user-images.githubusercontent.com/97848410/149844268-b6ca89b5-e124-490c-97c1-5163da8db24e.jpg)

Press no

11)	Wait till Well list refreshes (Will have a shortened Version of the original image name instead of A1-A12 and so on)

12)	Click on well and toggle on the Neurite Channel (corresponds to IBA1)

![Picture9](https://user-images.githubusercontent.com/97848410/149844269-e047311a-5dbf-48a6-9b64-4302481dadde.jpg)

13)	To lower intensity of the image, go to Options/Histogramm Options

14)	Set value to a higher number (e.g. 255) to lower brightness or set value to a lower number to make image brighter (Max Intensity Neurite for neurite channel and so on)

15)	Click on well name in well list to refresh

![Picture10](https://user-images.githubusercontent.com/97848410/149844270-61b2ad16-1447-4c51-8ed4-8e7f38bdd1a3.jpg)

16)	To draw regions, click on set filter and draw your polygon. Close the line by double clicking

17)	To delete regions, you need to go to the directory in which your raw images are saved and go to the Convertedcellomics folder. Within this folder you will find a matlab file which needs to be deleted (You just need to delete the file of the well for which you want to redraw the region!)

![Picture11](https://user-images.githubusercontent.com/97848410/149844271-719af416-8245-4737-8c64-612f58eeba8f.jpg)

18)	To run the analysis only on the ROI’s, go to Options/Neuronal Quantification Options and write a 1 into Use filters! (If you are not inserting a 1, the analysis will run very long, since the entire image will be analyzed)

![Picture12](https://user-images.githubusercontent.com/97848410/149844274-714fec4c-5f6f-4ef0-a76a-93a91c196265.jpg)

19)	To run the analysis navigate to Algorithms and select AnalyzeMicrogliaActivation

![Picture13](https://user-images.githubusercontent.com/97848410/149844275-5bfe2cbb-8aae-4490-a16f-b4f9676868c3.jpg)

20)	Enter an experiment name (this will be the final name of your excel file, so use something, which will help you to assign an excel file to an experiment!)

![Picture14](https://user-images.githubusercontent.com/97848410/149844276-547ae778-89ee-4bd4-98e7-f065364bfc15.jpg)

21)	Hit no for

![Picture15](https://user-images.githubusercontent.com/97848410/149844278-4d92e63b-be00-42fe-af53-0af6cef4e1e3.jpg)

22)	Hit no if you want to analyze your full ROI. If yes is chosen only a subpart of neurons I analzed and saved separatly for manual inspection.

![Picture16](https://user-images.githubusercontent.com/97848410/149844280-64100afe-b251-4f17-8f59-5c0011a742fe.jpg)

23)	You will now see a green progress bar. When it disappers, the analysis is completet and results are written into a excel file located in your Omnisphero driectoty with the name given in Experiment name!

24)	You can further save your analysis, by going to File/Save

![Picture17](https://user-images.githubusercontent.com/97848410/149844281-a92e86e7-c259-4436-b721-ba68ab74d773.jpg)

25)	Select a directory where you want to save the file and give it a name+.mat

![Picture18](https://user-images.githubusercontent.com/97848410/149844282-358d6a8e-37b7-4504-826a-592c08979f7b.jpg)

26)	The file will be saved and will appear as a zipped folder in your directory

![Picture19](https://user-images.githubusercontent.com/97848410/149844283-d036a07d-fab6-48b1-9baf-fc8e12337706.jpg)

27)	To reopen this file later, open Omnisphero as described above and go to File/Load

![Picture20](https://user-images.githubusercontent.com/97848410/149844284-74498852-e9ec-4083-9e71-60dda2dfe316.jpg)

28)	Select the .mat zip folder and press open (or double click on it)

29)	Wait till well list refreshes

30)	If you open the images from a path different under which you saved the images (e.g. if you open them on another computer) you need also to reload the image path. To do so, navigate to file 

![picture21](https://user-images.githubusercontent.com/97848410/149844285-c0b247f3-de98-4264-a72c-97d8261037ad.jpg)

31)	Select the parent image folder by just high lightening it and press select folder

![picture22](https://user-images.githubusercontent.com/97848410/149844286-e3347b9e-2e03-49b9-94df-7195e0f07c9c.jpg)

32)	Now you can select a well in list and toggle on and off channels!




Explanation sheets:

Sheet1:
Feature list on the single cell level (all cells, single and clustered cells)

Sheet2:
Histogram of ration total process length/cell body area

Sheet3:
Normalized Histogramm (Divided by number of fish)

Sheet4:
Summary for individual wells of total cell number, total area, the ROI area, the IBA1 Area IBA1 intensity, number of clustered cells

Sheet5:
Feature list on the single cell level (only single cells)

Sheet6:
Feature list on the single cell level (only clustered cells)

Sheet7:
Histogram of ration total process length/cell body area (only single cells)

Sheet8:
Histogram of ration total process length/cell body area (only clustered cells)


