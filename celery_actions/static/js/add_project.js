function alphanumeric(str) {
    if(str === "" || /^\s+$/.test(str)) {
        return "";
    }
    return str.match(/([0-9a-zA-Z ])/g).join("");
}

function toggleCloudButtonAction(){

    const cloudAddress = document.getElementById('cloudAddressInput');
    const toggleButton = document.getElementById('toggleCloudAddress');
    toggleButton.addEventListener('click', function () {
        if (cloudAddress.classList.contains('hidden')) {
            cloudAddress.classList.remove('hidden');
            toggleButton.innerText = 'Hide Cloud Address';
        } else {
            cloudAddress.classList.add('hidden');
            toggleButton.innerText = 'Show Cloud Address';
        }
    });
}

function addClassesFromText(){
    const separator=document.getElementById('separator').value;
    const allclasses=document.getElementById('classes').value;
    const classes = allclasses.split(separator);
    classes.forEach(function(c){
        if(c !== "" && /^\s+$/.test(c)===false){
            addClassItem(c);
        }
    })
    document.getElementById('classes').value=""

}

function addClassItem(name){

    function classItemInside(name){
        const valueInput = document.createElement('input');
        valueInput.className = 'project-class';
        valueInput.type = 'text';
        valueInput.value = name;
        valueInput.required = true;

        valueInput.addEventListener("input",function(e){
            e.target.value = alphanumeric(e.target.value);
        })

        const removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'btn btn-secondary removeExtra';
        removeButton.textContent = 'Remove';

        const extraGroup = document.createElement('div');
        extraGroup.appendChild(valueInput);
        extraGroup.appendChild(removeButton);

        return extraGroup;

    }
    const container = document.getElementById('classes_group');
    name=alphanumeric(name)
    if (name===""){
        return null;
    }
    const extraGroup = classItemInside(name)

    if (extraGroup!==null) {
        container.appendChild(extraGroup);
        extraGroup.querySelector('.removeExtra').addEventListener('click', function() {
            container.removeChild(extraGroup);
        });
    }else{
        console.log("there was an empty class");
    }


}

function addEmailsFromText(){
    const separator=";"
    const allemails=document.getElementById('emails').value;
    const emails = allemails.split(separator);
    emails.forEach(function(e){
        addEmailItem(e);
    })
    document.getElementById('emails').value=""
}

function addEmailItem(email){
    let inemail=email
    if(!email.includes("@")){
        inemail=email+"@nygenome.org"
    }
    function emailInside(email){
        const valueInput = document.createElement('input');
        valueInput.className = 'email';
        valueInput.type = 'email';
        valueInput.value = email;
        valueInput.required = true;

        const removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'btn btn-secondary removeExtra';
        removeButton.textContent = 'Remove';

        const extraGroup = document.createElement('div');
        extraGroup.appendChild(valueInput);
        extraGroup.appendChild(removeButton);

        return extraGroup;
    }
    const extraGroup = emailInside(inemail)
    const container = document.getElementById('emails_group');
    container.appendChild(extraGroup);
    extraGroup.querySelector('.removeExtra').addEventListener('click', function() {
        container.removeChild(extraGroup);
    });
}

function getProjectClasses(){
    const children=document.getElementById('classes_group').children;
    const classes=[]
    for (let i = 0; i < children.length; i++){
        const c=children[i];
        classes.push(c.children[0].value);
    }
    return classes;
}

function getProjectEmails(){
    const children=document.getElementById('emails_group').children;
    const emails=[]

    for (let i = 0; i < children.length; i++){
        const c=children[i];
        emails.push(c.children[0].value);
    }
    return emails;
}

function extraSectionHTML(name,type,value){
    function extraGroupInisde(name){
        // Create the input element for the name
        const nameInput = document.createElement('input');
        nameInput.type = 'text';
        nameInput.name = 'extraName';
        nameInput.placeholder = 'Name';
        nameInput.value = name
        nameInput.required = true;

        const typeSelect = document.createElement('select');
        typeSelect.name = 'extraType';
        typeSelect.required = true;

        const options = ['String', 'Boolean', 'Number'];
        options.forEach(optionText => {
            const option = document.createElement('option');
            option.value = optionText;
            option.textContent = optionText;
            typeSelect.appendChild(option);
        });
        typeSelect.value = type;

        const valueInput = document.createElement('input');
        valueInput.name = 'extraValue';

        if (typeSelect.value === 'Boolean') {
            valueInput.type = 'checkbox';
            valueInput.checked = false;
            if (typeof value === 'string' && value.includes('rue')) {
                valueInput.checked=true
            } else if (typeof value === 'boolean' && value === true) {
                valueInput.checked=true
            }
            valueInput.required=false;
        } else {
            valueInput.type = 'text';
            valueInput.value = ''; // Clear value when switching back to text
            valueInput.required=true;
        }
        valueInput.className = 'extra-group-value-input';
        valueInput.value = value
        valueInput.required = true;

        valueInput.addEventListener("input",function(e){
            if (e.target.type === 'text'){
                e.target.value = alphanumeric(e.target.value);
            }
        })

        const removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'btn btn-secondary removeExtra';
        removeButton.textContent = 'Remove';

        const container = document.createElement('div');
        container.appendChild(nameInput);
        container.appendChild(typeSelect);
        container.appendChild(valueInput);
        container.appendChild(removeButton);

        return container;
    }
    name= name || "";
    type= type || "String";
    value= value || "";

    const extraContainer = document.getElementById('extraContainer');
    const extraGroup=extraGroupInisde(name)
    extraGroup.classList.add('form-group', 'extra-group');
    extraContainer.appendChild(extraGroup);

    const extraTypeSelect = extraGroup.querySelector('select[name="extraType"]');
    const extraValueInput = extraGroup.querySelector('input[name="extraValue"]');

    function listenSetType(){
        if (extraTypeSelect.value === 'Boolean') {
            extraValueInput.type = 'checkbox';
            extraValueInput.checked = 'false'; // Default value for checkbox
            extraValueInput.required=false;
        } else {
            extraValueInput.type = 'text';
            extraValueInput.value = ''; // Clear value when switching back to text
            extraValueInput.required=true;
        }
    }

    extraTypeSelect.value=type
    listenSetType()
    extraValueInput.value = value

    extraTypeSelect.addEventListener('change',listenSetType);

    // Add event listener to remove button
    extraGroup.querySelector('.removeExtra').addEventListener('click', function() {
        extraContainer.removeChild(extraGroup);
    });

    return extraGroup;

}

function initExtraSection(){
    document.getElementById('addExtraBtn').addEventListener('click', function() {
        extraSectionHTML("","String","")
    });
}

function listenForSubmitFormAddProject(){
    const form = document.getElementById('dynamicForm');
    form.addEventListener('submit', function(event) {
        event.preventDefault();

        const emails=getProjectEmails();
        const classes=getProjectClasses();

        if(emails.length===0){
            alert('Please add emails');
            return;
        }
        if(classes.length===0){
            alert('Please add classes');
            return;
        }

        const formData = {
            cloud_address: document.getElementById('cloudAddressInput').value,
            token: document.getElementById('token').value,
            name: document.getElementById('name').value,
            rid: document.getElementById('rid').value,
            emails: emails,
            classes: classes,
            extra: []
        };

        document.querySelectorAll('.extra-group').forEach(group => {
            const extraName = group.querySelector('input[name="extraName"]').value;
            const extraType = group.querySelector('select[name="extraType"]').value;
            const extraValueInput = group.querySelector('input[name="extraValue"]');
            let extraValue;

            if (extraType === 'Boolean') {
                extraValue = extraValueInput.checked; // Get checkbox value
            } else {
                extraValue = extraValueInput.value; // Get text input value
            }

            formData.extra.push({ name: extraName, type: extraType, value: extraValue });
        });


        fetch('/add_project_do', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(formData)
        })
        .then(response => response.json())
        .then(data => {
            //alert(data.message);
            if ("error" in data) {
                document.getElementById('responseMessage').innerText = "Error within request"
                document.getElementById('responseData').innerText = data.error;
            }else{
                document.getElementById('responseMessage').innerText = "Sending this request"
                document.getElementById('responseData').innerText = JSON.stringify(data, null, 2);
            }

        })
        .catch(error => {
            console.log('Error:', error);
        });
    });

}

function initEventListeners(){
    document.getElementById('create_classes').addEventListener('click', addClassesFromText);
    document.getElementById('create_emails').addEventListener('click', addEmailsFromText);
}

function initView(){
    toggleCloudButtonAction();
    initExtraSection();
    extraSectionHTML("isVisium","Boolean",true);
    listenForSubmitFormAddProject();
    initEventListeners();
}

document.addEventListener('DOMContentLoaded', function () {
    initView();
});