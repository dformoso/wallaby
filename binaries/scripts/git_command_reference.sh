###########################################
############ FIRST TIME CONFIG ############
###########################################

git config --global user.name "daniel"
git config --global user.email "daniel.martinez.formoso@gmail.com"
git config --global color.ui true # Enable colored output in the terminal
git config --global core.editor emacs # Tells git that you want to use emacs
ssh-keygen -t rsa -C "daniel.martinez.formoso@gmail.com" # Set up ssh on your computer. Press Enter, Enter, Enter
cat ~/.ssh/id_rsa.pub # Copy your public key (the contents of the newly-created id_rsa.pub file) into your clipboard.
# go to: https://github.com/settings/keys - Add SSH Key” on the right. Add a label (like “My laptop”) and paste the public key into the big text box.
ssh -T git@github.com # Test your connection

############################################
################# wallaby #################
############################################

############### Do this once ###############
cd wallaby
git init
git remote add origin_remote git@github.com:dformoso/wallaby.git 
touch .gitignore

echo """
womtool-52.jar
cromwell-52.jar
diagrams
data/blast/
data/ref_genomes/
data/SRR5377828/
notebooks/.*
binaries/scripts/git_command_reference.txt
workflows/cromwell-executions/
workflows/cromwell-workflow-logs/
""" >> .gitignore

#############################################
####### Do this to commit new changes #######
#############################################

# There are three things to remember: "origin_remote" is one thing, "master" is another thing, and "origin_remote/master" is a separate thing. Three things total.
# - "master" is a local branch
# - "origin_remote/master" is a remote branch (which is a local copy of the branch named "master" on the remote named "origin_remote")
# "origin_remote" is a remote - not a branch

git branch -a # Display all and current branches
git remote show # Display all remotes
git remote show origin_remote  # Display info on specific remote
git status # Inspects the contents of the working directory and staging area
git add . # Adding everything to the staging area
git status # Inspects the contents of the working directory and staging area
git commit -m "changing readme"
git remote -v # Lists a Git project's remotes
git fetch origin_remote master # fetch master from the remote origin_remote. The master branch on origin will be fetched and the local copy will be named origin_remote/master
git merge origin_remote/master # Then merge origin_remote/master into master.
git push origin_remote master # Then you can push your new changes in master back to origin_remote
git log # Shows a list of all previous commits

#############################################
####### Do this to create a branch to #######
####### work on a new project feature #######
#############################################
###### INCOMPLETE ######
git checkout -b new_branch_name # Create the branch on your local machine and switch to this branch
git branch -a # Display all and current branches
git push origin_remote new_branch_name # Push the branch on github to remote called origin_remote
git branch -a # Display all and current branches
###### INCOMPLETE ######