/*
  ┬ ┬
  ├─┤TML and JS functionality
  ┴ ┴
*/

function make_element(type,obj){
    const e=document.createElement(type);
    for(k in obj){
        if(k.includes("class")){
            classes=obj[k].split(" ")
            for(c of classes){
                e.classList.add(c);
            }
        }else{
            e[k]=obj[k];
        }

    }
    return e;
}

function alphanumeric(str) {
    if(str === "" || /^\s+$/.test(str)) {
        return "";
    }
    return str.match(/([0-9a-zA-Z ])/g).join("");
}

/*
  ┌─┐
  │ ┬et DOM parts
  └─┘
*/

function getSRPath(){
    return document.getElementById("srdir").value
}

function getInput10XImageLocation(sampleid){
    const location=document.getElementById("10xloc_"+sampleid).value;
    return location;
}

function getInput10XImageScale(sampleid){
    let scale=document.getElementById("10xscale_"+sampleid).value;
    scale=Number.parseFloat(scale)
    return scale;
}

function getInputCustomImageLocation(sampleid){
    const location=document.getElementById("ownimgloc_"+sampleid).value;
    return location;
}

function getIdentifiedConvertedImageLocation(sampleid){
    const location=document.getElementById("customloc_"+sampleid).value;
    return location;
}

function getIdentifiedConvertedImageScale(sampleid){
    let scale=document.getElementById("customfinalscale_"+sampleid).value;
    if (scale===""){
        return ""
    }else{
        scale=Number.parseFloat(scale)
        return scale
    }
    return scale;
}

function getImageScaleBeforeConvert(sampleid){
    let scale=document.getElementById("ownimgscale_"+sampleid).value;
    if (scale===""){
        return ""
    }else{
        scale=Number.parseFloat(scale)
        return scale
    }

}

function getImageWidthBeforeConvert(sampleid){
    let width=document.getElementById("ownmax_"+sampleid).value;
    if (width===""){
        return ""
    }else{
        width=Number.parseInt(width)
        return width
    }
}

function getImageMaxWidthAutoBeforeConvert(sampleid){
    let width=document.getElementById("ownmaxauto_"+sampleid).value;
    if (width===""){
        return ""
    }else{
        width=Number.parseInt(width)
        return width
    }
}

function getImageCLI(sampleid){
    const box=document.querySelector("div[id=box_"+sampleid+"]")
    const cli=box.querySelector(".imgconsole");
    return cli;
}

function getSRCLI(sampleid) {
    const box = document.querySelector("div[id=box_" + sampleid + "]")
    const cli = box.querySelector(".srconsole");
    return cli;
}

function getDBCLI(sampleid){
    const box=document.querySelector("div[id=box_"+sampleid+"]")
    const cli=box.querySelector(".dbconsole");
    return cli;
}

function getGenesCLI(sampleid){
    const box=document.querySelector("div[id=box_"+sampleid+"]")
    const cli=box.querySelector(".geneconsole");
    return cli;
}

function getGeneConceptBoxsFilteredmatrixloc(sampleid){
    const box= document.querySelector("div[id=box_"+sampleid+"]")
    return box.querySelector(".filteredmatrixloc")
}

function getImageRadios(sampleid){
    const box=document.querySelector("div[id=box_"+sampleid+"]")
    const imageradios=box.querySelectorAll("[name=imagetype_"+sampleid+"]");
    return imageradios;
}

function clearSampleSelectionAreas(){
    const availableItems=document.getElementById("availableItems")
    availableItems.innerHTML=""
    const selectedItems=document.getElementById("selectedItems")
    selectedItems.innerHTML=""
    const matchimglocs=document.getElementById("matchimglocs")
    matchimglocs.innerHTML=""
    const samplesinfo=document.getElementById("samplesinfo")
    samplesinfo.innerHTML=""
}

function getSRfileDomInputs(sampleid){
    //{"selects":selects,"inputs":textpaths,"values":values,"paths":paths}
    const box=document.querySelector("div[id=box_"+sampleid+"]")
    const table=box.querySelector(".filestable");
    const children =table.children;
    let selects=[];
    let textpaths=[];
    let values=[];
    let paths=[];
    let checkups=[]
    for (let i=0;i<children.length;i++){
        if(i==0) continue;
        let row=children[i];
        selects.push(row.children[0].children[0]);
        textpaths.push(row.children[1].children[0]);
        values.push(row.children[0].children[0].value);
        paths.push(row.children[1].children[0].value);
        checkups.push(row.children[2].children[0].checked);
    }
    return {"selects":selects,"inputs":textpaths,"values":values,"paths":paths,"checks":checkups}
}

/*
  ┌─┐
  ├┤ etch functions
  └
*/

function fetchListSpacerangeDirs(){
    srdir=document.getElementById("srdir").value;
    if(srdir === "" || /^\s+$/.test(srdir)) {
        alert("Write a valid spaceranger directory")
    }

    clearSampleSelectionAreas();

    const formData = {
        spaceranger_dir: srdir
    };

    fetch('/list_sr_folders', {
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
                //react to error
            }else{

                availableItems.innerHTML=""

                const directories=data["directories"]
                for(const d of directories){
                    availableItems.innerHTML+=`<option value="${d}">${d}</option>`
                }
            }

        })
        .catch(error => {
            console.log('Error:', error);
        });



}

function fetchImageIdenitfy(sampleid,location){
    /*
    * sampleid: String
    * location: String ["10x"|"custom"]
    * */

    if(location!=="10x" && location!=="custom"){
        alert("Choose to identify 10X image or custom image")
        return;
    }

    let loc=""

    if(location==="10x"){
        loc=getInput10XImageLocation(sampleid);
    }else if(location==="custom"){
        loc=getInputCustomImageLocation(sampleid);
    }

    if(loc.length===0){
        alert("Fill the "+location+" image location to proceed.")
        return;
    }

    document.getElementById("uploadoption_"+sampleid).disabled=true

    const cli=getImageCLI(sampleid)
    cli.innerText="Asking server to identify "+location+" image properties ...\n"

    const formData = {"image_location": loc};
    fetch('/identify', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
        .then(response => response.json())
        .then(data => {
            astr=JSON.stringify(data)
            cli.innerText="success:\n"

            errorfound=false
            for(k of Object.keys(data)) {
                cli.innerText+=k+":"+data[k]+"\n"
                if(k.includes("error")){
                    errorfound=true
                }
            }

            if(errorfound===false){
                if (location==="10x"){
                    document.getElementById("option1_"+sampleid).checked=true
                    document.getElementById("option2_"+sampleid).checked=false
                    document.getElementById("uploadoption_"+sampleid).disabled=false
                }
            }

        })
        .catch(error => {
            cli.innerText="error locally:"+astr
        });
}

function fetchImageConvert(sampleid){

    // this is only called or custom image, get the scale,
    // maxwidth selected or if the max 9500 is selected

    //TODO: test more convert cases, check upload
    //test all the pipeline
    //make the fluorescence case,
    //make new repo, identify requirements

    const location=getInputCustomImageLocation(sampleid);
    const ownscalevalue= getImageScaleBeforeConvert(sampleid)
    const ownmaxwvalue = getImageWidthBeforeConvert(sampleid)
    const ownmaxwauto =getImageMaxWidthAutoBeforeConvert(sampleid)

    let width=ownmaxwauto
    let scale=1.0

    let sendWidth=true;
    let sendScale=false;

    if(ownscalevalue!=="" && ownmaxwvalue===""){
        sendScale=true
        sendWidth=false;
        scale=ownscalevalue
        console.log("send scale "+scale)
    }else if(ownmaxwvalue!=="" && ownscalevalue===""){
        let inwidth=ownmaxwvalue;
        if (inwidth<ownmaxwauto){
            width=inwidth
            document.getElementById("ownmax_"+sampleid).value=width
        }
        console.log("send width "+width)
    }

    let formData = {"image_location": location};
    if (sendWidth){
        formData["width"]=width
    }else if(sendScale){
        formData["scale"]=scale
    }

    const cli=getImageCLI(sampleid);

    cli.innerText="Asking server to convert image ...\n"

    //const formData = {"image_location": location,"image_type":imgclass};
    fetch('/convert_image', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
        .then(response => response.json())
        .then(data => {
            let cusloc=document.getElementById("customloc_"+sampleid);
            let cussca=document.getElementById("customfinalscale_"+sampleid);
            let upbutton=document.getElementById("uploadoption_"+sampleid);

            let errorfound=false

            astr=JSON.stringify(data)
            cli.innerText="success:\n"
            for(k of Object.keys(data)) {
                cli.innerText+=k+":"+data[k]+"\n"
                if (k.includes("error")){
                    errorfound=true
                }
            }

            if (cusloc)
                cusloc.value=data["convertloc"]
            if (cussca)
                cussca.value=data["scale"]
            if (upbutton)
                upbutton.disabled=false;

            if(errorfound===false){
                if (location==="10x"){
                    document.getElementById("option1_"+sampleid).checked=false
                    document.getElementById("option2_"+sampleid).checked=true
                    document.getElementById("uploadoption_"+sampleid).disabled=false
                }
            }

        })
        .catch(error => {
            cli.innerText="error locally:"+astr
        });
}

function fetchSearchFor10XHighResImage(sampleid){
    //cli
    const cli=getImageCLI(sampleid)
    cli.innerText="Searching for spaceranger related files...\n"
    //get scalefactors and srdir
    let srpath=getSRPath(sampleid);
    let files=getSRfileDomInputs(sampleid);
    let len=files["inputs"].length
    let scalefactorsfile=""
    for(let i=0;i<len;i++){
        if(files["paths"][i].includes("scalefactors_json.json")){
            scalefactorsfile=files["paths"][i]
        }
    }

    if(scalefactorsfile===""){
        alert("\"scalefactors_json.json\" file has not located, write the scale " +
            "manually or indicate the location of the \"scalefactors_json.json\" file")
    }

    const formData = {
        "sr_path":srpath,
        "sample":sampleid,
        "scalefactors_json":scalefactorsfile
    };

    const im10Xloc=document.getElementById("10xloc_"+sampleid);
    const im10Xscale=document.getElementById("10xscale_"+sampleid);

    fetch('/find_10x_highres_image', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
        .then(response => response.json())
        .then(data => {
            astr=JSON.stringify(data)
            cli.innerText="Response:\n"
            for(k of Object.keys(data)) {
                cli.innerText+=k+":"+data[k]+"\n"
                if(k.includes("image_loc") && data[k].length>0){
                    im10Xloc.value = data[k]
                }else if(k.includes("image_scale")){
                    im10Xscale.value = data[k]
                }
            }
        })
        .catch(error => {
            cli.innerText="error locally:"+astr
        });



}

function fetchSearchForSRFiles(sampleid){
    //cli
    const cli= getSRCLI(sampleid)
    cli.innerText="Searching for spaceranger related files...\n"
    //get list of items to search for, and also sr main path
    const srfilesinputs= getSRfileDomInputs(sampleid);
    const srpath=getSRPath();

    const formData = {
        "sr_path":srpath,
        "sample":sampleid,
        "file_list":srfilesinputs.values,
    };

    fetch('/find_srfiles', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
    .then(response => response.json())
    .then(data => {
        astr=JSON.stringify(data)
        cli.innerText="success:\n"
        for(k of Object.keys(data)) {
            cli.innerText+=k+":"+data[k]+"\n"
        }
        for(let i=0;i<data["file_locations"].length;i++){
            const filename=data["file_locations"][i][0]
            if(data["file_locations"][i][1].length>0){
                srfilesinputs["inputs"][i].value=data["file_locations"][i][1];

            }
            if(filename === "filtered_feature_bc_matrix.h5"){
                getGeneConceptBoxsFilteredmatrixloc(sampleid)
                    .value=data["file_locations"][i][1];
            }
            let len=data["file_locations"][i][1].length
            srfilesinputs["inputs"][i].scroll(len,0);
        }
        //console.log(data)

    })
    .catch(error => {
        cli.innerText="error locally:"+astr
    });
}

function fetchUploadImageToCloudRun(sampleid){
    //cli
    const cli=getImageCLI(sampleid)
    cli.innerText="Asking server to upload image to the cloud ...\n"
    //common for all
    const token = document.getElementById('token').value;
    if(token==="") {
        alert("Token is missing")
        cli.innerText="Google cloud token is missing ...\n"
        return;
    }
    const projectName = document.getElementById('projectName').value;
    if(projectName==="") {
        alert("Project name is missing")
        cli.innerText="Project name is missing.\n"
        return;
    }

    const rid = document.getElementById('rid').value;
    if(rid==="") {
        alert("Role is missing")
        cli.innerText="Role is missing.\n"
        return;
    }

    let option="option1"
    if(document.getElementById("option2_"+sampleid).checked){
        option="option2"
    }

    let isIF =document.getElementById("isifccheck_"+sampleid).checked
    let group=null;
    let marker=null;

    if(isIF){
        group=document.getElementById('ifcmarkergroup_'+sampleid).value;
        marker=document.getElementById('ifcmarker_'+sampleid).value;
    }

    let imagePath=""
    let scale=""

    if (option.includes("option1")){
        //specific to sample
        imagePath = getInput10XImageLocation(sampleid);
        scale = getInput10XImageScale(sampleid);
    }else if(option.includes("option2")){
        imagePath = getIdentifiedConvertedImageLocation(sampleid);
        scale = getIdentifiedConvertedImageScale(sampleid);
    }



    const formData = {
        "token":token,
        "image_path":imagePath,
        "rid":rid,
        "sample":sampleid,
        "project":projectName,
        "scale":scale
    };

    if (isIF){
        formData["is_IF"]=true;
        formData["IF_group"]=group;
        formData["IF_marker"]=marker;
    }

    fetch('/upload_image', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
        .then(response => response.json())
        .then(data => {
            //alert(data.message);
            astr=JSON.stringify(data)
            cli.innerText="success:\n"
            for(k of Object.keys(data)) {
                cli.innerText+=k+":"+data[k]+"\n"
            }
        })
        .catch(error => {
            console.log('Error:', error);
        });
}

function fetchUploadGenesForSample(sampleid){
    //cli
    const cli=getGenesCLI(sampleid)
    cli.innerText="Asking server to upload selected genes ...\n"

    const token = document.getElementById('token').value;
    if(token==="") {
        alert("Token is missing")
        return;
    }
    const projectName = document.getElementById('projectName').value;
    const rid = document.getElementById('rid').value;

    //get the genes from the box for that
    const geneboxvalue=document.getElementById('genelist').value;

    if(geneboxvalue.length===0){
        alert("No genes selected");
        return
    }

    let genes=geneboxvalue.split(/[,|\n]/);


    let trimmedlist = genes.filter(str => str.trim() !== "");

    console.log("genes: "+trimmedlist)

    //get the is ensembl
    const isensembl=document.getElementById("isensembl_"+sampleid).checked;

    const srdominputs=getSRfileDomInputs(sampleid)

    console.log(srdominputs);

    let tissuepositions=""
    let filteredmatrixh5=""
    for(let i=0;i<srdominputs["values"].length;i++){
        let val=srdominputs["values"][i] || ""
        let path=srdominputs["paths"][i] || ""
        if(val.includes("tissue") && val.includes(".csv")){
            tissuepositions=path
        }
        if(val.includes("filtered") && val.includes(".h5")){
            filteredmatrixh5=path
        }
    }

    console.log("tissuepositions ,filteredmatrixh5,",tissuepositions ,filteredmatrixh5)

    if(tissuepositions.length ===0 || filteredmatrixh5.length ===0){
        alert("The necessary files were not found")
        return
    }
        //{"selects":selects,"inputs":textpaths,"values":values,"paths":paths}

    // get the tissue positions and the filtered matrix h5

    const formData = {
        "token":token,
        "rid":rid,
        "sample":sampleid,
        "project":projectName,
        "genes":genes,
        "tissue_positions_loc":tissuepositions,
        "filtered_matrix_h5_loc":filteredmatrixh5,
        "is_ensembl":isensembl
    };

    fetch('/upload_gene_data', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
        .then(response => response.json())
        .then(data => {
            //alert(data.message);
            astr=JSON.stringify(data)
            cli.innerText="success:\n"
            for(k of Object.keys(data)) {
                let val=data[k]+"\n"
                cli.innerText+=k+":"+val
            }
        })
        .catch(error => {
            console.log('Error:', error);
        });

}

function fetchUploadSpots(sampleid){

    //get locations of relevant files
    const srfiles=getSRfileDomInputs(sampleid);

    tp_path="";metadata_path=""; ms_path="";

    for(i=0;i<srfiles["values"].length;i++){
        let v=srfiles["values"][i];
        if(v.includes("tissue_positions.csv")){
            tp_path=srfiles["paths"][i];
        }
        if(v.includes("scalefactors_json.json")){
            metadata_path=srfiles["paths"][i];
        }
        if(v.includes("metrics_summary.csv")){
            ms_path=srfiles["paths"][i];
        }
    }

    if(tp_path==="") {
        alert("tissue positions path must exist")
        return;
    }

    const token = document.getElementById('token').value;
    if(token==="") {
        alert("Token is missing")
        return;
    }
    const projectName = document.getElementById('projectName').value;
    const rid = document.getElementById('rid').value;

    //specific to sample
    const sampleName = sampleid;

    const formData = {
        "token":token,
        "rid":rid,
        "sample":sampleName,
        "project":projectName,
        "tp_path":tp_path
    };

    if(metadata_path !== "") {
        formData["metadata_path"] = metadata_path
    }
    if(metadata_path !== "") {
        formData["ms_path"] = ms_path
    }

    const cli=getSRCLI(sampleid)
    cli.innerText="Asking server to upload sample and spot infos ...\n"

    fetch('/upload_spot_info', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
        .then(response => response.json())
        .then(data => {
            astr=JSON.stringify(data)
            cli.innerText="success:\n"
            for(k of Object.keys(data)) {
                cli.innerText+=k+":"+data[k]+"\n"
            }
        })
        .catch(error => {
            cli.innerText="error locally:"+astr
        });
}

function fetchCreateSampleInDB(sampleid){
    //cli
    const cli=getDBCLI(sampleid)
    cli.innerText="Asking server to insert in database ...\n"
    //common for all
    const token = document.getElementById('token').value;
    if(token==="") {
        alert("Token is missing")
        return;
    }
    const projectName = document.getElementById('projectName').value;
    const rid = document.getElementById('rid').value;

    //specific to sample
    const sampleName = sampleid;

    const formData = {
        "token":token,
        "rid":rid,
        "sample":sampleName,
        "project":projectName,
    };

    fetch('/initialize_in_db', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
    .then(response => response.json())
    .then(data => {
        //alert(data.message);
        astr=JSON.stringify(data)
        cli.innerText="success:\n"
        for(k of Object.keys(data)) {
            let val=data[k]+"\n"
            // this is to avoid confusion as no xy barcode is uploaded only an empty text file,
            // we just need to create the item in the database , this is because of how the cloud run function
            // is programmed, we could change it there too, but that would be for someone else to handle, it is
            // not critical
            val=val.replace("and uploaded xy_barcode json file to storage bucket","")
            cli.innerText+=k+":"+val
        }
    })
    .catch(error => {
        console.log('Error:', error);
    });

}

/*
  ┌─┐
  │  reate interface parts
  └─┘
*/

function createDBConceptBox(sampleid){
    const dbboxconcept=make_element("div",{"class":"sampleboxconcept"});
    const dbc_title=make_element("div",{"innerText":"Create in database:","class":"title"});
    const dbc_explain=make_element("p", {"innerText": "If the sample doesn't exist in the database it has " +
                                                          "to be created first. Make sure the project name and the role id " +
                                                      "are correct and that the project exists!"})//,{"class":""});
    const dbc_facconsole=make_element("p",{"class":"facconsole dbconsole","innerText":""});

    const dbc_createindb=make_element("button",{"type":"button","id":"createindbbtn_"+sampleid,"class":"createindbbtn","innerText":"Create entry in database"});

    dbc_createindb.addEventListener("click",(event)=>{
        const sampleid=event.currentTarget.id.replace("createindbbtn_","");
        fetchCreateSampleInDB(sampleid);
    });

    dbboxconcept.appendChild(dbc_title);
    dbboxconcept.appendChild(dbc_explain);
    dbboxconcept.appendChild(dbc_createindb);
    dbboxconcept.appendChild(dbc_facconsole);

    return dbboxconcept;
}

function _toggleMarkerHide(event){
    event.preventDefault();
    let target=event.target;
    const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
    const elem=document.getElementById("markerififcbox_"+parent_id)
    if (target.checked) {
        elem.classList.remove("hidden")
    }else{
        elem.classList.add("hidden")
    }
}

function createImageConceptBox(sampleid){
    const imageboxconcept=make_element("div",{"class":"sampleboxconcept"});

    const imc_title=make_element("div",{"innerText":"Image data:","class":"title"});

    //---- option 1 --------------------------------------
    const imc_option1box=make_element("div",{"class":"imoptcontain"});
    const imc_liloptcolumnmain=make_element("div",{"class":"imoptcolumn"});
    const imc_labopt1=make_element("span",{"class":"imoptrottext","innerText":"Option 1"});
    const imc_radioopt1=make_element("input",{"type":"radio" ,"id":"option1_"+sampleid,"name":"optioninput_"+sampleid, "value":"option1", "class":"radio"})
    const imc_restoption1container=make_element("div",{"class":"imoptrestoption"});
    const imc_choosehighrestitle=make_element("div",{"innerText":"10x's tissue_hires_image.png"});
    const imc_10xhighresloc=make_element("input",{"type":"text", "style":"width:inherit","id":"10xloc_"+sampleid});
    const imc_10xscaletitle=make_element("div",{"innerText":"Scale according to scale factors json " +
                                                            "(if available. Make sure the file is selected in the" +
                                                            "space ranger files box above otherwise copy manually)"});
    const imc_10xhighreslscale=make_element("input",{"type":"text", "style":"width:35%","id":"10xscale_"+sampleid});
    const imc_buttonrow=make_element("div",{"class":"imoptbuttonrow"});
    const imc_button10xfind=make_element("button",{"type":"button","id":"10xfind_"+sampleid,"class":"findimagebtn","innerText":"find"});
    const imc_button10xidentify=make_element("button",{"type":"button","id":"10xidentify_"+sampleid,"class":"identifyimagebtn","innerText":"identify"});

    imc_button10xfind.addEventListener("click",event=>{
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        fetchSearchFor10XHighResImage(parent_id);
    });

    imc_button10xidentify.addEventListener("click",event=>{
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        fetchImageIdenitfy(parent_id,"10x");
    })

    imc_radioopt1.checked=true;

    imc_liloptcolumnmain.appendChild(imc_radioopt1);
    imc_liloptcolumnmain.appendChild(imc_labopt1);
    imc_buttonrow.appendChild(imc_button10xfind);
    imc_buttonrow.appendChild(imc_button10xidentify);
    imc_restoption1container.appendChild(imc_choosehighrestitle);
    imc_restoption1container.appendChild(imc_10xhighresloc);
    imc_restoption1container.appendChild(imc_10xscaletitle);
    imc_restoption1container.appendChild(imc_10xhighreslscale);
    imc_restoption1container.appendChild(imc_buttonrow);
    imc_option1box.appendChild(imc_liloptcolumnmain);
    imc_option1box.appendChild(imc_restoption1container);

    //---- option 2 --------------------------------------

    const imc_option2box=make_element("div",{"class":"imoptcontain"});
    const imc_liloptcolumnmain2=make_element("div",{"class":"imoptcolumn"});
    const imc_labopt2=make_element("span",{"class":"imoptrottext","innerText":"Option 2"});
    const imc_radioopt2=make_element("input",{"type":"radio" ,"id":"option2_"+sampleid,"name":"optioninput_"+sampleid, "value":"option2", "class":"radio"})
    const imc_restoption2container=make_element("div",{"class":"imoptrestoption"});

    const imc_ownimgtitle=make_element("div",{"innerText":"Location of custom image"});
    const imc_ownimgloc=make_element("input",{"type":"text", "style":"width:inherit" ,"id":"ownimgloc_"+sampleid});

    const imc_ownscalerow=make_element("div",{"class":"scalerowcontainer","style":"width:inherit"});
    const imc_ownscalecol=make_element("div",{"class":"colscale","innerHTML":"<span>Scale</span>"});
    const imc_ownscalecoltext=make_element("input",{"type":"text", "style":"","id":"ownimgscale_"+sampleid});
    //const imc_ownor1=make_element("div",{"style":"margin:20px 10px 0 10px;","innerHTML":"or"});
    const imc_ownor1=make_element("div",{"style":"margin-top:17px","innerHTML":"<span>or</span>"});
    const imc_ownmaxwcol=make_element("div",{"class":"colscale","innerHTML":"<span>Max width (max 9500)</span>"});
    const imc_ownmaxwcoltext=make_element("input",{"type":"text", "style":"","id":"ownmax_"+sampleid});
    const imc_ownor2=make_element("div",{"style":"margin-top:17px","innerHTML":"<span>or</span>"});
    const imc_ownautomaxwcol=make_element("div",{"class":"colscale","innerHTML":"<span>Auto (max default width)</span>"});
    const imc_ownautomaxwcoltext=make_element("input",{"type":"text", "value":"9500","id":"ownmaxauto_"+sampleid});
    imc_ownautomaxwcoltext.disabled=true;
    const imc_buttonrow2=make_element("div",{"class":"imoptbuttonrow"});
    const imc_buttonownimgidentify=make_element("button",{"type":"button","id":"ownidentify_"+sampleid,"class":"identifyimagebtn","innerText":"identify"});

    const imc_createcustomrow=make_element("div",{"class":"scalerowcontainer","style":"width:inherit"});
    const imc_customfinalscalelabel=make_element("div",{"class":"colscale","innerHTML":"<span>Final scale</span>"});
    const imc_customfinalscale=make_element("input",{"type":"text", "style":"","id":"customfinalscale_"+sampleid});
    const imc_customloclabel=make_element("div",{"class":"colscale","style":"flex:3","innerHTML":"<span>Location of image to upload</span>"});
    const imc_customloc=make_element("input",{"type":"text", "style":"","id":"customloc_"+sampleid});

    const imc_buttoncreateim=make_element("button",{"type":"button","id":"owncreate_"+sampleid,"class":"createownimbtn","innerText":"Create with these settings"});

    imc_buttoncreateim.addEventListener("click",event=>{
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        fetchImageConvert(parent_id,"custom");
    })

    imc_buttonownimgidentify.addEventListener("click",event=>{
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        fetchImageIdenitfy(parent_id,"custom");
    })

    imc_radioopt1.addEventListener("change",event=>{
            const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
            document.getElementById("uploadoption_"+parent_id).disabled=true
    });
    imc_radioopt2.addEventListener("change",event=>{
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        document.getElementById("uploadoption_"+parent_id).disabled=true
    })

    imc_liloptcolumnmain2.appendChild(imc_radioopt2);
    imc_liloptcolumnmain2.appendChild(imc_labopt2);
    imc_ownscalecol.appendChild(imc_ownscalecoltext);
    imc_ownmaxwcol.appendChild(imc_ownmaxwcoltext);
    imc_ownautomaxwcol.appendChild(imc_ownautomaxwcoltext);
    imc_ownscalerow.appendChild(imc_ownscalecol);
    imc_ownscalerow.appendChild(imc_ownor1);
    imc_ownscalerow.appendChild(imc_ownmaxwcol);
    imc_ownscalerow.appendChild(imc_ownor2);
    imc_ownscalerow.appendChild(imc_ownautomaxwcol);
    imc_buttonrow2.appendChild(imc_buttonownimgidentify);
    imc_buttonrow2.appendChild(imc_buttoncreateim);
    imc_customfinalscalelabel.appendChild(imc_customfinalscale)
    imc_customloclabel.appendChild(imc_customloc);
    imc_createcustomrow.appendChild(imc_customfinalscalelabel)
    imc_createcustomrow.appendChild(imc_customloclabel)
    //imc_createbuttonrow.appendChild(imc_buttoncreateim);
    imc_restoption2container.appendChild(imc_ownimgtitle);
    imc_restoption2container.appendChild(imc_ownimgloc);
    imc_restoption2container.appendChild(imc_ownscalerow);
    imc_restoption2container.appendChild(imc_buttonrow2);
    imc_restoption2container.appendChild(imc_createcustomrow);
    //imc_restoption2container.appendChild(imc_createbuttonrow);
    imc_option2box.appendChild(imc_liloptcolumnmain2);
    imc_option2box.appendChild(imc_restoption2container);

    const imc_isvisifcrow=make_element("div",{"class":"rowcontainer"});
    const imc_isvisiumcheck=make_element("input",{"type":"checkbox","id":"isvisiumcheck_"+sampleid});
    const imc_isvisiumchecklabel=make_element("label",{"for":"isvisiumcheck_"+sampleid, "innerText":"is visium","class":"blacklabel"});
    const imc_isifccheck=make_element("input",{"type":"checkbox","id":"isifccheck_"+sampleid});
    const imc_isifcchecklabel=make_element("label",{"for":"isifccheck_"+sampleid, "innerText":"is fluorescence","class":"blacklabel"});

    imc_isvisifcrow.appendChild(imc_isifccheck);
    imc_isvisifcrow.appendChild(imc_isifcchecklabel);
    imc_isvisifcrow.appendChild(imc_isvisiumcheck);
    imc_isvisifcrow.appendChild(imc_isvisiumchecklabel);

    const imc_markerififc=make_element("div",{"id":"markerififcbox_"+sampleid,"class":"hidden"});
    const imc_markergrouptitle=make_element("div",{"innerText":"Group"});
    const imc_markergroupttext=make_element("input",{"type":"text", "style":"width:35%","id":"ifcmarkergroup_"+sampleid});
    const imc_markertitle=make_element("div",{"innerText":"Marker"});
    const imc_markertext=make_element("input",{"type":"text", "style":"width:35%","id":"ifcmarker_"+sampleid});

    imc_markerififc.appendChild(imc_markergrouptitle);
    imc_markerififc.appendChild(imc_markergroupttext);
    imc_markerififc.appendChild(imc_markertitle);
    imc_markerififc.appendChild(imc_markertext);

    const imc_facconsole=make_element("p",{"class":"facconsole imgconsole","innerText":""});

    const imc_uploadbtnrow=make_element("div",{"class":"rowcontainer"});
    const imc_buttoncuploadselected=make_element("button",{"type":"button","id":"uploadoption_"+sampleid,"innerText":"Upload"});
    imc_buttoncuploadselected.disabled=true;
    imc_uploadbtnrow.appendChild(imc_buttoncuploadselected);

    imc_buttoncuploadselected.addEventListener("click",event=>{
        event.preventDefault();
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","");
        fetchUploadImageToCloudRun(parent_id);
    })

    imc_isifccheck.addEventListener("input", _toggleMarkerHide)
    imc_isifccheck.addEventListener("change", _toggleMarkerHide)

    imageboxconcept.appendChild(imc_title);
    imageboxconcept.appendChild(imc_option1box);
    imageboxconcept.appendChild(imc_option2box);
    imageboxconcept.appendChild(imc_isvisifcrow);
    imageboxconcept.appendChild(imc_markerififc);
    imageboxconcept.appendChild(imc_uploadbtnrow);
    imageboxconcept.appendChild(imc_facconsole);

    return imageboxconcept;

}

function createSpacerangerFileInput(option){
    const src_row=make_element("tr");
    const src_td1=make_element("td");
    const src_td2=make_element("td");
    const src_td3=make_element("td");
    const src_td4=make_element("td");

    const src_tdselect=make_element("select",{"class":"srfile"});
    let options=["barcodes.tsv.gz", "matrix.mtx.gz", "features.tsv.gz", "tissue_positions.csv",
        "tissue_positions_list.csv","tissue_positions.parquet",
        "scalefactors_json.json","metrics_summary.csv","filtered_feature_bc_matrix.h5","other"];
    options.forEach(op => {
        const opt = make_element("option",{"value":op,"innerText":op});
        if (op === option) opt.selected = true;
        src_tdselect.appendChild(opt);
    })
    src_td1.appendChild(src_tdselect);
    const src_tdinputtext=make_element("input",{"type":"text","class":"srlocation"});
    src_td2.appendChild(src_tdinputtext);
    const src_tdinputcheck=make_element("input",{"type":"checkbox","class":"srupload"});
    src_td3.appendChild(src_tdinputcheck);
    const src_tdremovebtn=make_element("button",{"type":"button","class":"removerow","innerText":"X"});
    src_td4.appendChild(src_tdremovebtn);

    src_tdremovebtn.addEventListener("click", (event)=>{
        event.preventDefault();
        const target=event.target;
        const rowparent=target.parentElement.parentElement; //get row, not td
        rowparent.parentElement.removeChild(rowparent);
    });

    src_row.appendChild(src_td1)
    src_row.appendChild(src_td2)
    src_row.appendChild(src_td3)
    src_row.appendChild(src_td4)

    return src_row;

}

function createSpacerangerConceptBox(sampleid){
    const srcboxconcept=make_element("div",{"class":"sampleboxconcept"});
    const src_title=make_element("div",{"innerText":"Spaceranger data:","class":"title"});
    const src_table=make_element("table",{"class":"filestable"});
    const src_tableheader=make_element("tr",{"innerHTML":"<th>sr file</th>  <th>location</th> <th>upload</th> <th>remove</th>"});

    src_table.appendChild(src_tableheader);

    let options=["tissue_positions.csv", "filtered_feature_bc_matrix.h5", "scalefactors_json.json", "metrics_summary.csv", "barcodes.tsv.gz", "matrix.mtx.gz", "features.tsv.gz",  "other"];
    options.forEach(op => {
        const row = createSpacerangerFileInput(op);
        src_table.appendChild(row);
    });

    const src_rowbtncontainer=make_element("div",{"class":"rowcontainer"})

    const src_autofindbtn=make_element("button",{"type":"button","id":"srautofind_"+sampleid,"class":"srautofind","innerText":"Auto find all"});
    src_rowbtncontainer.appendChild(src_autofindbtn);
    const src_addfilebtn=make_element("button",{"type":"button","class":"addsrbtn","innerText":"Add SR file"});
    src_rowbtncontainer.appendChild(src_addfilebtn);
    src_addfilebtn.addEventListener("click",(event)=>{
        const row=createSpacerangerFileInput("other");
        src_table.appendChild(row);
    });
    const src_upallbtn=make_element("button",{"type":"button","style":"margin-left:auto;","class":"sruploadall","innerText":"Mark upload all"});
    src_rowbtncontainer.appendChild(src_upallbtn);

    srcboxconcept.appendChild(src_title);
    srcboxconcept.appendChild(src_table);
    srcboxconcept.appendChild(src_rowbtncontainer);

    //new container only for the "create spots" if visium
    const src_rowcreatecontainer=make_element("div",{"class":"rowcontainer","style":"margin-top:5px;"})
    const src_ifvp=make_element("span",{"innerText":"If Visium: ", "style":"padding:0px 10px 0px 10px"});
    const src_createspotsifvbtn=make_element("button",{"type":"button","id":"createspots_"+sampleid,"class":"createspots","innerText":"Create spots"});
    const src_needp=make_element("span",{"innerText":"(needs the location of 'tissue_positions.csv') ", "style":"padding:0px 5px 0px 7px; color:#aaaaaa"});
    src_rowcreatecontainer.appendChild(src_ifvp);
    src_rowcreatecontainer.appendChild(src_createspotsifvbtn);
    src_rowcreatecontainer.appendChild(src_needp);
    srcboxconcept.appendChild(src_rowcreatecontainer);

    src_createspotsifvbtn.addEventListener("click",(event)=>{
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        fetchUploadSpots(parent_id);
    });

    const src_facconsole=make_element("p",{"class":"facconsole srconsole","innerText":""});
    srcboxconcept.appendChild(src_facconsole);

    //const src_upsrdatabtn=make_element("button",{"type":"button","class":"uploadsr","innerText":" Upload Spaceranger data"});
    //srcboxconcept.appendChild(src_upsrdatabtn);

    src_upallbtn.addEventListener("click",(event)=>{
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        const upchecks = src_table.querySelectorAll(".srupload");
        if(src_upallbtn.innerText.includes("Unmark")){
            upchecks.forEach((ch)=>{
                ch.checked=false;
            })
            src_upallbtn.innerText="Mark upload all"
        }else{
            upchecks.forEach((ch)=>{
                ch.checked=true;
            })
            src_upallbtn.innerText="Unmark upload all"
        }

    });

    src_autofindbtn.addEventListener("click",(event)=>{
        //get all the items in the box and find them
        //const sampleid=event.currentTarget.id.replace("srautofind_","");
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        fetchSearchForSRFiles(parent_id)

    });

    return srcboxconcept;
}

function createGeneConceptBox(sampleid){
    const geneboxconcept=make_element("div",{"class":"sampleboxconcept"});
    //<label for="vehicle1"> I have a bike</label><br>
    const gcb_title=make_element("div",{"class":"title","innerText":"Genes for visualization:"});
    const gcb_part1=make_element("div",{"innerText":"Gene list: Will be the same as the general gene list" +
                                                      " selected in the previous sections"})
    const gcb_isensemblexplain=make_element("div",{"innerText":"This marks if the genes have ensembl names"});
    const gcb_isensembllabel=make_element(("label"),{"for":"isensembl_"+sampleid, "innerText":"is ensembl",
    "style":"color:black; font-weight: normal;"});
    const gcb_isensembl=make_element("input",{"id":"isensembl_"+sampleid,"type":"checkbox","class":"gcisensembl"});
    gcb_isensembl.checked=false;
    const gcb_part2=make_element("div",{"innerText":"H5 gene filtered matrix location. Same as the one " +
                                                        "found in the space ranger files in the box before this one"})
    const gcb_text=make_element("input",{"type":"text","size":80,"class":"filteredmatrixloc"})
    const gcb_facconsole=make_element("p",{"class":"facconsole geneconsole","innerText":""});
    const gcb_uploadgenes=make_element("button",{"type":"button","class":"uploadparquet","innerText":"Upload gene data","style":"display: block;"});

    geneboxconcept.appendChild(gcb_title);
    geneboxconcept.appendChild(gcb_part1);
    geneboxconcept.appendChild(gcb_isensemblexplain);
    geneboxconcept.appendChild(gcb_isensembl);
    geneboxconcept.appendChild(gcb_isensembllabel);
    geneboxconcept.appendChild(gcb_part2);
    geneboxconcept.appendChild(gcb_text);
    geneboxconcept.appendChild(gcb_facconsole);
    geneboxconcept.appendChild(gcb_uploadgenes);

    gcb_uploadgenes.addEventListener("click",(event)=>{
        //get all the items in the box and find them
        //const sampleid=event.currentTarget.id.replace("srautofind_","");
        const parent_id=event.target.closest(".inputsamplebox").id.replace("box_","")
        fetchUploadGenesForSample(parent_id)

    });

    return geneboxconcept;
}

function createInputRegion(sampleid) {

    const inputRegions = document.getElementById('samplesinfo');

    const inputsamplebox=make_element('div',{"id":"box_"+sampleid,"class":"inputsamplebox"});

    const sampleboxlabel=make_element ('label',{"innerText":sampleid});
    const collapsea=make_element("a",{"class":"collapse","innerText":">","href":"#"})

    const wholesamplecontainer=make_element('div',{"class":"wholesamplecontainer hidden"});
    // ^ contains all concept boxes

    sampleboxlabel.appendChild(collapsea);
    inputsamplebox.appendChild(sampleboxlabel);
    inputsamplebox.appendChild(wholesamplecontainer);

    collapsea.addEventListener("click", (event)=>{
        event.preventDefault();
        const target = event.target;
        const parent = event.target.parentElement;
        if (wholesamplecontainer.classList.contains("hidden")) {
            wholesamplecontainer.classList.remove("hidden")
            target.innerText="^";
        }else{
            wholesamplecontainer.classList.add("hidden")
            target.innerText=">";
        }
    });

    const dbboxconcept=createDBConceptBox(sampleid);
    wholesamplecontainer.appendChild(dbboxconcept);

    const srboxconcept=createSpacerangerConceptBox(sampleid);
    wholesamplecontainer.appendChild(srboxconcept);

    const imageboxconcept=createImageConceptBox(sampleid);
    wholesamplecontainer.appendChild(imageboxconcept);

    const geneboxconvept=createGeneConceptBox(sampleid);
    wholesamplecontainer.appendChild(geneboxconvept);

    inputsamplebox.appendChild(wholesamplecontainer);
    inputRegions.appendChild(inputsamplebox);
}

function reactToIFC(target){
    const name=target.getAttribute("name");
    const sampleid=target.closest(".inputsamplebox").id.replace("box_","")
    const imageradios=document.querySelectorAll("[name="+name+"]")
    const markercontainer = document.getElementById("markercontainer_"+sampleid)
    console.log(markercontainer)
    imageradios.forEach(imageradio=>{
        if(imageradio.checked){
            if(imageradio.classList.contains("wsiifc") || imageradio.classList.contains("visiumifc")){
                markercontainer.classList.remove("hidden");
            }else{
                markercontainer.classList.add("hidden");
            }
        }
    });
}

function addContentBoxesForSelected(){
    const availableItems = document.getElementById('availableItems');
    const selectedItems = document.getElementById('selectedItems');
    const inputRegions = document.getElementById('inputRegions');

    // Move selected items from availableItems to selectedItems
    for (let i = 0; i < availableItems.options.length; i++) {
        if (availableItems.options[i].selected) {
            const option = availableItems.options[i];
            selectedItems.add(new Option(option.text, option.value));
            availableItems.remove(i);
            i--; // Decrement i to account for the removed item

            // Create a new input region for the selected item
            createInputRegion(option.value);
        }
    }
}

function removeContentBoxesFromSelected(){
    const selectedItems = document.getElementById('selectedItems');
    const availableItems = document.getElementById('availableItems');
    const inputRegions = document.getElementById('samplesinfo');

    // Move selected items from selectedItems back to availableItems
    for (let i = 0; i < selectedItems.options.length; i++) {
        if (selectedItems.options[i].selected) {
            const option = selectedItems.options[i];
            availableItems.add(new Option(option.text, option.value));
            selectedItems.remove(i);
            i--; // Decrement i to account for the removed item

            // Remove the corresponding input region
            const inputRegion = document.getElementById("box_"+option.value);
            if (inputRegion) {
                inputRegions.removeChild(inputRegion);
            }
        }
    }
}

function copyImageLocationsToIndividualBoxes(){
    const locationbox=document.getElementById("matchimglocs");
    const selectedbox=document.getElementById("selectedItems");

    const imlocs=locationbox.value.split("\n")

    if(imlocs.length !== selectedbox.options.length){
        alert("Amount of selected samples and corresponding image locations must match")
        return;
    }

    for(let i=0; i<selectedbox.options.length;i++){
        const sampname=selectedbox.options[i];
        const location=imlocs[i]
        if(location !== ""){
            const box=document.getElementById("box_"+sampname);
            const imtext=box.querySelector('.pathtoimage');
            imtext.value=location;
        }
    }
}

/*
  ┬
  │  isten, external functions
  ┴─┘
*/

function listenSpacerangeDirsPress(){
    document.getElementById("listsrdir").addEventListener("click", fetchListSpacerangeDirs)
    document.getElementById("selectedAddBtn").addEventListener("click", addContentBoxesForSelected)
    document.getElementById("selectedRemoveBtn").addEventListener("click", removeContentBoxesFromSelected)

}

function listenCopyLocationsToSamples(){
    document.getElementById("copyimagelocsbtn").addEventListener("click", copyImageLocationsToIndividualBoxes);

}

function listenSamplesInfoButtons(){
    document.getElementById("allaswsibtn").addEventListener("click", ()=>{
        const wsiradios=document.querySelectorAll('input.wsibright.radio');
        wsiradios.forEach(radio => {
            radio.checked=true;
        })
    });

    document.getElementById("allasvisiumbtn").addEventListener("click", ()=>{
        const visiumradios=document.querySelectorAll('input.visiumbright.radio');
        visiumradios.forEach(radio => {
            radio.checked=true;
        })
    })

    document.getElementById("allasifbtn").addEventListener("click", ()=>{
        const ifradios=document.querySelectorAll('input.wsiifc.radio');
        ifradios.forEach(radio => {
            radio.checked=true;
        });
        ifradios.forEach(radio => {
            reactToIFC(radio)
        });
    })

    document.getElementById("allasifbtn").addEventListener("click", ()=>{
        const ifradios=document.querySelectorAll('input.visiumifc.radio');
        ifradios.forEach(radio => {
            radio.checked=true;
        });
        ifradios.forEach(radio => {
            reactToIFC(radio)
        });
    })

    document.getElementById("collapseall").addEventListener("click", (event)=>{
        wholesamplecontainers=document.querySelectorAll(".wholesamplecontainer");
        if(event.target.innerText.includes("Collapse")){
            wholesamplecontainers.forEach(elem=>{
                elem.classList.add("hidden")
            });
            event.target.innerText="Show all"
        }else{
            wholesamplecontainers.forEach(elem=>{
                elem.classList.remove("hidden")
            })
            event.target.innerText="Collapse all"
        }
    })

}

function fillDummy() {
    document.getElementById("srdir").value = "/home/leslie";//"/commons/groups/vickovic_lab/crc/visium/spaceranger_outputs/2024-07-22_run/"
    //document.getElementById("prefixpathimage").value = "/gpfs/commons/groups/vickovic_lab/gridnext/data/human_ba46/images/ccast_removed/V087_V10S15-032-A1.jpg"
    document.getElementById("projectName").value="exampleproject";
    document.getElementById("rid").value="CGND";

}

function initView(){

    fillDummy();

    listenSpacerangeDirsPress();
    listenForSubmitFormAddProject();
    listenSamplesInfoButtons();
    listenCopyLocationsToSamples();
}

function listenForSubmitFormAddProject(){
    const form = document.getElementById('dynamicForm');
    form.addEventListener('submit', function(event) {
        event.preventDefault();
    });
}

document.addEventListener('DOMContentLoaded', initView );