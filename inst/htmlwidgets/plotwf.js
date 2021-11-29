HTMLWidgets.widget({
  name: 'plotwf',
  type: 'output',
  factory: function(el, width, height) {
    // TODO: define shared variables for this instance
    /*functions */
    function makeResponsive(plot){
     if (plot) {
      if (plot.width) {plot.removeAttribute("width")}
      if (plot.height) {plot.removeAttribute("height")}
      plot.style.width = "100%";
      plot.style.height = "100%";
     }
    }
    function toPng(){
      let plot_el = el.querySelector('svg');
      domtoimage.toBlob(plot_el)
        .then(function (blob) {
            window.saveAs(blob, 'plotwf.png');
      });
    }

    function toJpg(){
      let plot_el = el.querySelector('svg');
        domtoimage.toJpeg(plot_el, { quality: 1 })
    .then(function (dataUrl) {
        var link = document.createElement('a');
        link.download = 'plotwf.jpeg';
        link.href = dataUrl;
        link.click();
    });
    }

    function toSvg(){
      let plot_el = el.querySelector('svg');
      domtoimage.toSvg(plot_el)
      .then(function (dataUrl) {
        var link = document.createElement('a');
        link.download = 'plotwf.svg';
        link.href = dataUrl;
        link.click();
      });
    }
    //function toSvg(){
    //  let plot_el = el.querySelector('svg');
    //  var file = new File([plot_el.outerHTML], "plotwf.svg", {type: "text/plain;charset=utf-8"});
    //  saveAs(file);
    //}

    function toPdf(){
      let plot_el = el.querySelector('svg');
      domtoimage.toPng(plot_el)
      .then(function(dataurl) {
          var orientation = window.innerHeight >= window.innerWidth ? "portrait" : "landscape";
          const doc = new jspdf.jsPDF({unit: "px", orientation: orientation, format: [window.innerHeight, window.innerWidth], compress: true});
          doc.addImage(dataurl, "PNG", 0, 0, window.innerWidth, window.innerHeight);
          doc.save("plotwf.pdf");
      });
    }

    function toDot(dotStr){
      var file = new File([dotStr], "plotwf.dot", {type: "text/plain;charset=utf-8"});
      saveAs(file);
    }

    async function load_scripts(script_urls) {
        function load(script_url) {
            return new Promise(function(resolve, reject) {
                if (load_scripts.loaded.has(script_url)) {
                    resolve();
                } else {
                    var script = document.createElement('script');
                    script.onload = resolve;
                    script.src = script_url
                    document.head.appendChild(script);
                }
            });
        }
        var promises = [];
        for (const script_url of script_urls) {
            promises.push(load(script_url));
        }
        await Promise.all(promises);
        for (const script_url of script_urls) {
            load_scripts.loaded.add(script_url);
        }
    }
    load_scripts.loaded = new Set();

    function addControl(el, dotStr) {
      (async () => {
      await load_scripts([
        'https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.4.0/jspdf.umd.min.js',
        'https://cdnjs.cloudflare.com/ajax/libs/FileSaver.js/2.0.0/FileSaver.min.js'
        ]);

        var ctrGroup = document.createElement("div")
        ctrGroup.className = "wfplot-ctr"

        if(!domtoimage || !jspdf || !saveAs) {
          ctrGroup.innerHTML = `<p class="error-msg">JS libraries not loaded, check your Internet</p>`;
          el.appendChild(ctrGroup);
          return false;
        } else {

        }


        ctrGroup.innerHTML =
        `
        <button data-for="png">PNG</button>
        <button data-for="jpg">JPG</button>
        <button data-for="svg">SVG</button>
        <button data-for="pdf">PDF</button>
        <button data-for="dot">Graphiz</button>
        `
        el.appendChild(ctrGroup);
        ctrGroup.querySelector('button[data-for="png"]').addEventListener('click', toPng)
        ctrGroup.querySelector('button[data-for="jpg"]').addEventListener('click', toJpg)
        ctrGroup.querySelector('button[data-for="svg"]').addEventListener('click', toSvg)
        ctrGroup.querySelector('button[data-for="pdf"]').addEventListener('click', toPdf)
        ctrGroup.querySelector('button[data-for="dot"]').addEventListener('click', ()=>{toDot(dotStr)})
      })();

      return el
    }
    /**************/

    return {
      renderValue: function(x) {
        // don't know why Rmd like to append this to table and doesn't give it a </div> closing, remove it
        if (x.rmd) x.dot = x.dot.replaceAll(`<div class='horizontal-scroll'>`, '');
        var plot = document.querySelector(`#${x.plotid}`);

        var viz = new Viz();
        var legendSrc = document.querySelector('head link[id*="plotwf_legend"]').attributes.href.value;
        var dotStr = x.dot.replace('plotwf_legend-src\.png', legendSrc);
        [...el.children].forEach(e => e.remove()); //clear all children before attach

        viz[x.plot_method](dotStr, {images: [{path: legendSrc, width: '450px', height: '250px'}]})
        .then(function(plot_el) {
          plot_el.id = x.plotid;
          el.style.width = x.width ? x.width: "100%";
          el.style.height = x.height ? x.height: "100%";
          el.style.overflow = "auto";

          if(x.msg !== "") {
            let msgNode = document.createElement("h4");
            msgNode.innerText = x.msg;
            msgNode.style.textAlign = "center";
            msgNode.style.color = "red";
            el.appendChild(msgNode);
          }
          el.appendChild(plot_el);
          var singlePage = document.querySelector(`#htmlwidget_container`);
          if(singlePage) document.querySelector('body').style.padding = '0';

          if (x.rmd) $(plot_el).removeAttr('height width');

          if(x.responsive) makeResponsive(plot_el);
          document.dispatchEvent(new Event('wf_plot_created'));

          return el;
        })
        .catch(e => {
          var p = document.createElement("pre");
          console.log(e);
          p.innerText = e;
          el.appendChild(p);
        });

        if(x.plot_ctr) addControl(el, dotStr);

      },
      resize: function(width, height) {
        // TODO: code to re-render the widget with a new size
      }
    };
  }
});


