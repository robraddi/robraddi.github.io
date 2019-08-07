---
layout: post
title:  "git-branch & git-merge"
author: Rob
categories: [ Git ]
featured: true
---


# Quick Note:
### **Always commit with a useful message.** Otherwise, your history will be useless...

```bash
# Show your most recent commits
git show-branch --all
```

Assuming prior knowledge of
```bash
git pull; git add; git commit; git push
```

### \#<span style="color:red">TODO</span>:

- **Fork** a repository: [octocat/Spoon-Knife](https://github.com/octocat/Spoon-Knife)

- Create a new **branch**

- Make changes to existing code or develop code.

- **Merge** your new branch to the master branch to create a **Pull Request**


# After forking octocat/Spoon-Knife — get a local copy:

```bash
git clone https://github.com/USERNAME/Spoon-Knife.git
```

## Creating a new branch

<img src="https://wac-cdn.atlassian.com/dam/jcr:fcad863b-e0da-4a55-92ee-7caf4988e34e/02.svg?cdnVersion=le" width="400" height="400" align="right">


```bash
git branch [new branch]
git checkout [new branch]      # git checkout -b [new branch]
git add [stuff/work]
git commit -m "my new work"
git push --set-upstream origin [new branch]
```

### Check what's been done:

```bash
git show-branch
```

## Merging new branch to another branch
```bash
git checkout [branch]
git merge [new branch]
git push origin [branch]
```

### Now, in your browser create a **pull request**.
<br>

<footer> https://www.atlassian.com/git/tutorials/why-git</footer>


## Fetching the Pull Request (For Maintainer):

## Step 1: From your project repository, bring in the changes and test.
```bash
git fetch origin
git checkout -b [new branch] origin/[new branch]
git merge master
```
## Step 2: Merge the changes and update on GitHub.
```bash
git checkout master
git merge --no-ff [new branch]
git push origin master
```



# Other Notes:


### Information commands to learn about git

```bash
man giteveryday

man gittutorial
```



