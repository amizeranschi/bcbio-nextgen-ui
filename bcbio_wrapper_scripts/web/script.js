const defaultPages = ["run_config", "downstream_report", "multiqc_report", "fastqc_report", "about", "help"]


/* -----------------------------------------------------------------FUNCTIONS FOR RUN CONFIGURATION PAGE ---------------------------------------------------------*/

// Used the create the appearance of the multi-page application
// This function will interchange and load html pages in index.html file
// The navbar buttons will trigger this function and load one of the pages from web folder.
function loadPage(page, style) {

	// First disable all the pages from the available list
	defaultPages.forEach(pageTemplate => {
		document.getElementById(pageTemplate).style["display"] = "none";
	})

	// Than load the required one
	document.getElementById(page).style["display"] = "block";
	document.getElementById(page).innerHTML = '<object type="text/html" style="' + style + '" data="./' + page + '.html" ></object>'
}

function showSpinner() {
	parent.document.getElementById("spinner").style["display"] = "block";
}

function hideSpinner() {
	parent.document.getElementById("spinner").style["display"] = "none";
}

// Used to hide on show some input fields based on a value from checkbox
/*
	For example: existing_version - yes 
					hide - install_path field
					show - path_to_existing
*/
function displayFieldOnCheckBox(id, onTrueFields, onFalseFields) {
	// Get the checkbox
	var checkBox = document.getElementById(id);

	// We can have multiple fields to show/display
	// So we will iterate through all of them and show/hide based on the checkbox value
	onTrueFields.forEach(field => {
		document.getElementById(field).style.display = checkBox.checked == true ? "block" : "none"
	})
	onFalseFields.forEach(field => {
		document.getElementById(field).style.display = checkBox.checked == false ? "block" : "none"
	})
}

// Populate in-line error under the input fields
function printError(elemId, hintMsg) {
	const err_field = document.getElementById(elemId + '_err')
	if (err_field)
		err_field.innerHTML = hintMsg;
}

// Used to display pop-ups on the page
/* For this function we have 2 types of modals
	- confirmation
	- error
*/
function openModal(elemId, type, title, message, style) {

	// Get the modal - all the modals will be in the same page, but they will be hidden
	var modal = document.getElementById(elemId);

	// Get the <span> element that closes the modal
	var closeButtons = document.getElementsByClassName("close-span");

	// Display the modal and set the other properties
	modal.style.display = "block";
	document.getElementById("modal-msg_" + type).innerHTML = message;
	document.getElementById("modal-title_" + type).innerHTML = title;

	// Add the style
	const modalContents = document.querySelectorAll('.modal-content');
	modalContents.forEach(modalContent => {
		modalContent.style.cssText = style;
	});

	// When the user clicks on <span> (x), close the modal
	for (let index = 0; index < closeButtons.length; index++) {
		closeButtons[index].onclick = function () {
			modal.style.display = "none";
		}
	}

	// When the user clicks anywhere outside of the modal, close it
	window.onclick = function (event) {
		if (event.target == modal) {
			modal.style.display = "none";
		}
	}
}

// Used to display the Preview Modal of the YAML file
function openPreviewModal(content, style) {
	// Get the modal
	var modal = document.getElementById("preview-modal");
	// Get the <span> element that closes the modal
	var closeButtons = document.getElementsByClassName("close-span");

	// Display the modal and set the other properties
	modal.style.display = "block";
	document.getElementById("preview-yml").innerHTML = content;

	// Add the style
	const modalContents = document.querySelectorAll('.modal-content-top');
	modalContents.forEach(modalContent => {
		modalContent.style.cssText = style;
	});

	// When the user clicks on <span> (x), close the modal
	for (let index = 0; index < closeButtons.length; index++) {
		closeButtons[index].onclick = function () {
			modal.style.display = "none";
		}
	}

	// When the user clicks anywhere outside of the modal, close it
	window.onclick = function (event) {
		if (event.target == modal) {
			modal.style.display = "none";
		}
	}
}

// Generate the HTML content based on the YAML file and trigger the show modal function
function generatePreview() {
	eel.generate_preview()(function (content) {
		openPreviewModal(content, 'width: 70vw')
	})
}

// Validate the Form and check if all the required fields are populated and display error message in case they are not
function validateForm(element) {
	hasErrors = false;
	// We need to get the parent div of the field to check if is on hidden mode or not.
	// Because there is no need to validate a filed if is hidden
	var divOfElement = document.getElementById(element.id + '_div')

	if (!divOfElement || (divOfElement && (divOfElement.style.display == 'block' || divOfElement.style.display == ''))) {
		// Different validation flows based on the element type

		if (element.type == 'checkbox' && (element.checked == null || element.checked == undefined)) {
			printError(element.id, 'Please fill out this field')
			hasErrors = true;
		} else if (element.type !== 'checkbox' && (element.value == null || element.value == undefined || element.value == "")) {
			printError(element.id, 'Please fill out this field')
			hasErrors = true;
		} else {
			printError(element.id, '')
		}
	}

	return hasErrors
}

// Load an existing configuration file from local system
function loadFile() {
	hasErrorsList = []

	// Open File System Browse modal and select a file
	eel.load_config_file()(function (content) {
		// Based on the result from python function we will either populate the flow or throw an error
		if (content) {
			// No action needs to be taken if the user can cancel the File System Browse modal 
			if (content != 'cancel') {
				const form = document.getElementById('form');
				const formElements = Array.from(form.elements);
				formElements.forEach(element => {
					if (!element.type.startsWith('submit') && !element.type.startsWith('button')) {
						if (element.type == 'checkbox') {
							element.checked = content[element.id] == 'yes' ? true : false;
							if (typeof element.onclick == "function") {
								element.onclick()
							}
						} else {
							element.value = content[element.id]
						}

						// We still want to validate the flow in case loaded file is incomplete
						hasErrors = validateForm(element)
						if (hasErrors)
							hasErrorsList.push(hasErrors)
					}
				})

				document.getElementById('preview-btn').disabled = !(hasErrorsList.length == 0);
				document.getElementById('run-analysis-btn').disabled = !(hasErrorsList.length == 0);

				if (hasErrorsList.length > 0) {
					openModal('error-modal', 'err', 'Errors Found in Configuration Form',
						'<p>There are fields that require your attention.<p/> <p>Please fill out all required fields!<p/>', 'width: 40%')
				}
			}

		} else {
			// Throw an error in case user is trying to load a wrong file type
			openModal('error-modal', 'err', 'Invalid File',
				'<p>Please check the extension and the content of the file.<p/> <p>Only .yml and .yaml extension files are accepted.<p/>', 'width: 50%')
		}
	})
}


// Used to display reports dropdown menu only after the user ran analysis
function runAnalysis() {
	parent.document.getElementsByClassName("dropdown")[0].style.visibility = "visible"
	parent.document.getElementsByClassName("tip-below")[0].style.display = "none"

	showSpinner()
	eel.run_analysis()(function (result) {
		const workflow_elem = document.getElementById("workflow")
		hideSpinner()
		if (result) {
			openModal('confirmation-modal', 'succ', 'Successfully ran the ' + (workflow_elem ? workflow_elem.value : "") + " Workflow",
			'<p>You can check the the results on the Report Tab.<p/>', 'width: 50%')
		} else {
			openModal('error-modal', 'err', 'The ' + (workflow_elem ? workflow_elem.value : "") + "Workflow Failed",
				'<p>Lorem Ipsum Error message.<p/>', 'width: 50%')
		}
	})
}

function handleSubmit(e) {
	e.preventDefault();

	var object = {}
	var hasErrorsList = [];
	const form = document.getElementById('form');

	// Get all form elements and look for ids
	const formElements = Array.from(form.elements);

	formElements.forEach(element => {

		// Validate Form
		if (!element.type.startsWith('submit') && !element.type.startsWith('button')) {
			hasErrors = validateForm(element)
			if (hasErrors)
				hasErrorsList.push(hasErrors)
		}

		// Populate JSON
		if (element.type.startsWith('select')) {
			const re = /^\d*(\.\d+)?$/
			if (element.value.match(re)) {
				object[element.id] = { type: 'number', value: element.value };
			} else {
				object[element.id] = { type: 'text', value: element.value };
			}
		} else if (element.type == 'checkbox') {
			object[element.id] = { type: element.type, value: element.checked ? "yes" : "no" };
		} else if (!element.type.startsWith('submit') && !element.type.startsWith('button')) {
			object[element.id] = { type: element.type, value: element.value };
		}
	});

	if (hasErrorsList.length == 0) {
		eel.save_configuration(object)(function (result) {
			if (result != '') {
				openModal('confirmation-modal', 'succ', 'Successfully saved YAML file',
					'<p>The configuration file was saved at the following location.<p/> <p style="word-break: break-all">' + result + '<p/>', 'width: 50%')
				document.getElementById('preview-btn').disabled = false;
				document.getElementById('run-analysis-btn').disabled = false;
			}
			else
				openModal('error-modal', 'err', 'Unable to save YAML file',
					'<p>There are fields that require your attention.<p/> <p>Please fill out all required fields!<p/>', 'width: 40%')
		})
	}
	else
		openModal('error-modal', 'err', 'Errors Found in Configuration Form',
			'<p>There are fields that require your attention.<p/> <p>Please fill out all required fields!<p/>', 'width: 40%')
}

const form = document.querySelector('form');
if (form)
	form.addEventListener('submit', handleSubmit);


/* -----------------------------------------------------------------FUNCTIONS FOR DOWNSTREAM REPORTS ---------------------------------------------------------*/

function generateHtmlComponents(components) {

	components_st = ''
	components.forEach(component => {
		components_st +=
				'<div class="component-report">' +
					'<div>' +
						'<span style="white-space: pre-wrap;">'
							+ component["description"] +
						'</span>' +
					'</div>' +
					'<div>' +
						'<img style="width: 100%;" src="' + component["imagePath"] + '" />' +
					'</div>' +
				'</div>'
	})

	return components_st;
}
function generateHtmlReports() {
	const main_container = document.getElementById('main-container');
	var html_st = ""
	eel.get_report_context_file()(function (reportContext) {
		// Add Workflow Name
		html_st += '<h1 class="workflow-title">' + reportContext['workflow'] + '</h1>'
	
		// Add content
		reportContext['content'].forEach(content => {
			// Add Sub title
			html_st += 
					'<div class="container-box">' + 
						'<h3 class="sub-title">' + content["title"] + '</h3>' +
						generateHtmlComponents(content["components"]) +
					'</div>'
		})
	
		main_container.innerHTML = html_st
	})
}

