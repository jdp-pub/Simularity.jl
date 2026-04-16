# Introduction

Work in progress (WIP). 

Thank you for your interest in Simularity, the high performance computing simulation toolkit for Julia. External developers give the toolkit its refined edge, as it is iron sharpens iron. The purpose of this document is to create a procedure as the project grows in size and scope while also preserving maintainability for the forseeable future.

## Accepted Contributions

Pretty much anthing relevant will be accepted as long as its in the appropriate form. Bug fixes, new functionality pertaining to simulations, optomizations, code generalizations, documentation fixes are all wanted. What is not wanted is anything that is broken or breaks something else (outside of the working branches). Partially complete functionality or the beginning of optomizations or anything similar is welcome in working branches.


# Workflow

##Branches:

`stable` - current release

`staging` - next release, stable code

`<feature,bugfix,optomizations,another category>\<topic>_working` - unstable workspace, make this when nonexistent. whenever recent changes are complete sync with staging branch. Some features may take longer before staging is appropriate.

## Checklist before staging

- unit tests
- functions are arranged in lexicographical order per file
- [complete code documentation](https://docs.julialang.org/en/v1/manual/documentation/)
- update documentation pages (TBA)

[julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/#Write-functions,-not-just-scripts)

[performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/)

## Conventions

This project uses snake_case for names in general, variable names may also be standard mathematical notation like A for matrices or x for vectors/scalars, i,j,k, for indeces, etc. as long as the purpose is apparent based on documentation and function.

## Dependencies

The project is not "native only" but dependencies will be carefully considered for necessity. If dependency overhead can reasonably be reduced within the scope of the toolkit then it will be rewritten and incorporated. 

## Depth

Please limit function depth. Nothing should be repeated but unnecessary function depth to make code more "readable" is actually a form of obfuscation that hinders performance. 

## Size

Additions can be any size, but working branches should be complete before staging.

# Documentation

## Changelog format

Include a hyperlink to function documentation any time a function is mentioned.



[format directions](https://keepachangelog.com/en/1.1.0/)

## Citations

Cite all claims, anyhting original should be published and then cited. APS formatting, hyperlinked to source.

#### Formats

ArXiv: `[^key]: A. Author, arXiv:(arxiv number) (Year).`
Book: `[^key]: A. Author, *Book Title* (Publisher, City, Year).`
Conference/Proceedings: `[^key]: A. Author, in *Proceedings Title* (Publisher, City, Year), pp. (start page)-(end page).`
Journal Article: `[^key]: A. Author, B. Author, *Journal Abbrev.* **Volume**, Page (Year).`
Thesis: `[^key]: A. Author, Ph.D. thesis, University Name (Year).`
Websites: `[^key]: Title or Organization, URL (accessed Month Day, Year).`


<!--
template items, will be incorporated if theres any interest in the project

### Bonus points: Add a link to a resource for people who have never contributed to open source before.
Here are a couple of friendly tutorials you can include: http://makeapullrequest.com/ and http://www.firsttimersonly.com/

> Working on your first Pull Request? You can learn how from this *free* series, [How to Contribute to an Open Source Project on GitHub](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github).

[source: [React](https://github.com/facebook/react/blob/master/CONTRIBUTING.md#pull-requests)]  

As a side note, it helps to use newcomer-friendly language throughout the rest of your document. Here are a couple of examples from [Active Admin](https://github.com/activeadmin/activeadmin/blob/master/CONTRIBUTING.md):

>At this point, you're ready to make your changes! Feel free to ask for help; everyone is a beginner at first :smile_cat:
>
>If a maintainer asks you to "rebase" your PR, they're saying that a lot of code has changed, and that you need to update your branch so it's easier to merge.

# Getting started
### Give them a quick walkthrough of how to submit a contribution.
How you write this is up to you, but some things you may want to include:

* Let them know if they need to sign a CLA, agree to a DCO, or get any other legal stuff out of the way
* If tests are required for contributions, let them know, and explain how to run the tests
* If you use anything other than GitHub to manage issues (ex. JIRA or Trac), let them know which tools they’ll need to contribute

>For something that is bigger than a one or two line fix:

>1. Create your own fork of the code
>2. Do the changes in your fork
>3. If you like the change and think the project could use it:
    * Be sure you have followed the code style for the project.
    * Sign the Contributor License Agreement, CLA, with the jQuery Foundation.
    * Note the jQuery Foundation Code of Conduct.
    * Send a pull request indicating that you have a CLA on file.

[source: [Requirejs](http://requirejs.org/docs/contributing.html)] **Need more inspiration?** [1] [Active Admin](https://github.com/activeadmin/activeadmin/blob/master/CONTRIBUTING.md#1-where-do-i-go-from-here) [2] [Node.js](https://github.com/nodejs/node/blob/master/CONTRIBUTING.md#code-contributions) [3] [Ember.js](https://github.com/emberjs/ember.js/blob/master/CONTRIBUTING.md#pull-requests)

### If you have a different process for small or "obvious" fixes, let them know.

> Small contributions such as fixing spelling errors, where the content is small enough to not be considered intellectual property, can be submitted by a contributor as a patch, without a CLA.
>
>As a rule of thumb, changes are obvious fixes if they do not introduce any new functionality or creative thinking. As long as the change does not affect functionality, some likely examples include the following:
>* Spelling / grammar fixes
>* Typo correction, white space and formatting changes
>* Comment clean up
>* Bug fixes that change default return values or error codes stored in constants
>* Adding logging messages or debugging output
>* Changes to ‘metadata’ files like Gemfile, .gitignore, build scripts, etc.
>* Moving source files from one directory or package to another

[source: [Chef](https://github.com/chef/chef/blob/master/CONTRIBUTING.md#chef-obvious-fix-policy)] **Need more inspiration?** [1] [Puppet](https://github.com/puppetlabs/puppet/blob/master/CONTRIBUTING.md#making-trivial-changes)

# How to report a bug
### Explain security disclosures first!
At bare minimum, include this sentence:
> If you find a security vulnerability, do NOT open an issue. Email XXXX instead.

If you don’t want to use your personal contact information, set up a “security@” email address. Larger projects might have more formal processes for disclosing security, including encrypted communication. (Disclosure: I am not a security expert.)

> Any security issues should be submitted directly to security@travis-ci.org
> In order to determine whether you are dealing with a security issue, ask yourself these two questions:
> * Can I access something that's not mine, or something I shouldn't have access to?
> * Can I disable something for other people?
>
> If the answer to either of those two questions are "yes", then you're probably dealing with a security issue. Note that even if you answer "no" to both questions, you may still be dealing with a security issue, so if you're unsure, just email us at security@travis-ci.org.

[source: [Travis CI](https://github.com/travis-ci/travis-ci/blob/master/CONTRIBUTING.md)] **Need more inspiration?** [1] [Celery](https://github.com/celery/celery/blob/master/CONTRIBUTING.rst#security) [2] [Express.js](https://github.com/expressjs/express/blob/master/Security.md)

### Tell your contributors how to file a bug report.
You can even include a template so people can just copy-paste (again, less work for you).

> When filing an issue, make sure to answer these five questions:
>
> 1. What version of Go are you using (go version)?
> 2. What operating system and processor architecture are you using?
> 3. What did you do?
> 4. What did you expect to see?
> 5. What did you see instead?
> General questions should go to the golang-nuts mailing list instead of the issue tracker. The gophers there will answer or ask you to file an issue if you've tripped over a bug.

[source: [Go](https://github.com/golang/go/blob/master/CONTRIBUTING.md#filing-issues)] **Need more inspiration?** [1] [Celery](https://github.com/celery/celery/blob/master/CONTRIBUTING.rst#other-bugs ) [2] [Atom](https://github.com/atom/atom/blob/master/CONTRIBUTING.md#reporting-bugs) (includes template)

# How to suggest a feature or enhancement
### If you have a particular roadmap, goals, or philosophy for development, share it here.
This information will give contributors context before they make suggestions that may not align with the project’s needs.

> The Express philosophy is to provide small, robust tooling for HTTP servers, making it a great solution for single page applications, web sites, hybrids, or public HTTP APIs.
>
> Express does not force you to use any specific ORM or template engine. With support for over 14 template engines via Consolidate.js, you can quickly craft your perfect framework.

[source: [Express](https://github.com/expressjs/express#philosophy)] **Need more inspiration?** [Active Admin](https://github.com/activeadmin/activeadmin#goals)

### Explain your desired process for suggesting a feature.
If there is back-and-forth or signoff required, say so. Ask them to scope the feature, thinking through why it’s needed and how it might work.

> If you find yourself wishing for a feature that doesn't exist in Elasticsearch, you are probably not alone. There are bound to be others out there with similar needs. Many of the features that Elasticsearch has today have been added because our users saw the need. Open an issue on our issues list on GitHub which describes the feature you would like to see, why you need it, and how it should work.

[source: [Elasticsearch](https://github.com/elastic/elasticsearch/blob/master/CONTRIBUTING.md#feature-requests)] **Need more inspiration?** [1] [Hoodie](https://github.com/hoodiehq/hoodie/blob/master/CONTRIBUTING.md#feature-requests) [2] [Ember.js](https://github.com/emberjs/ember.js/blob/master/CONTRIBUTING.md#requesting-a-feature)

# Code review process
### Explain how a contribution gets accepted after it’s been submitted.
Who reviews it? Who needs to sign off before it’s accepted? When should a contributor expect to hear from you? How can contributors get commit access, if at all?

> The core team looks at Pull Requests on a regular basis in a weekly triage meeting that we hold in a public Google Hangout. The hangout is announced in the weekly status updates that are sent to the puppet-dev list. Notes are posted to the Puppet Community community-triage repo and include a link to a YouTube recording of the hangout.
> After feedback has been given we expect responses within two weeks. After two weeks we may close the pull request if it isn't showing any activity.

[source: [Puppet](https://github.com/puppetlabs/puppet/blob/master/CONTRIBUTING.md#submitting-changes)] **Need more inspiration?** [1] [Meteor](https://meteor.hackpad.com/Responding-to-GitHub-Issues-SKE2u3tkSiH ) [2] [Express.js](https://github.com/expressjs/express/blob/master/Contributing.md#becoming-a-committer)

# Community
If there are other channels you use besides GitHub to discuss contributions, mention them here. You can also list the author, maintainers, and/or contributors here, or set expectations for response time.

> You can chat with the core team on https://gitter.im/cucumber/cucumber. We try to have office hours on Fridays.

[source: [cucumber-ruby](https://github.com/cucumber/cucumber-ruby/blob/master/CONTRIBUTING.md#talking-with-other-devs)] **Need more inspiration?**
 [1] [Chef](https://github.com/chef/chef/blob/master/CONTRIBUTING.md#-developer-office-hours) [2] [Cookiecutter](https://github.com/audreyr/cookiecutter#community)

# BONUS: Code, commit message and labeling conventions
These sections are not necessary, but can help streamline the contributions you receive.

### Explain your preferred style for code, if you have any.

**Need inspiration?** [1] [Requirejs](http://requirejs.org/docs/contributing.html#codestyle) [2] [Elasticsearch](https://github.com/elastic/elasticsearch/blob/master/CONTRIBUTING.md#contributing-to-the-elasticsearch-codebase)

### Explain if you use any commit message conventions.

**Need inspiration?** [1] [Angular](https://github.com/angular/material/blob/master/.github/CONTRIBUTING.md#submit) [2] [Node.js](https://github.com/nodejs/node/blob/master/CONTRIBUTING.md#step-3-commit)

### Explain if you use any labeling conventions for issues.

**Need inspiration?** [1] [StandardIssueLabels](https://github.com/wagenet/StandardIssueLabels#standardissuelabels) [2] [Atom](https://github.com/atom/atom/blob/master/CONTRIBUTING.md#issue-and-pull-request-labels)

-->
