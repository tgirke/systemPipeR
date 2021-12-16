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

    function toSvgFull(){
      let plot_el = el.querySelector('svg');
      domtoimage.toSvg(plot_el)
      .then(function (dataUrl) {
        var link = document.createElement('a');
        link.download = 'plotwf_full_size.svg';
        link.href = dataUrl;
        link.click();
      });
    }

    function toSvgMin(){
      let plot_el = el.querySelector('svg');
      var file = new File([plot_el.outerHTML], "plotwf_min_size.svg", {type: "text/plain;charset=utf-8"});
      saveAs(file);
    }

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

    async function load_scripts(script_urls, scripts_loaded) {
        function load(script_url) {
            return new Promise(function(resolve, reject) {
                if (scripts_loaded.has(script_url)) {
                    resolve();
                } else {
                    var script = document.createElement('script');
                    script.onload = resolve;
                    script.src = script_url;
                    script.addEventListener('load', ()=>scripts_loaded.add(script_url))
                    document.head.appendChild(script);
                }
            });
        }
        var promises = [];
        for (const script_url of script_urls) {
            promises.push(load(script_url));
        }
        await Promise.allSettled(promises);
    }

    function ctrInernetErr(ctrGroup, el) {
      let errMsg = "Some JS libraries not loaded, check your Internet";
      ctrGroup.innerHTML = `<p class='error-msg'>${errMsg}</p>`;
      el.appendChild(ctrGroup);
      throw new Error(errMsg);
    }

    function panZoomInernetErr(panErr, el) {
      let errMsg = "Pan-zoom is not loaded, check your Internet";
      panErr.innerHTML = `<p class='error-msg'>${errMsg}</p>`;
      el.appendChild(panErr);
      throw new Error(errMsg);
    }

    function addControl(el, dotStr) {
      (async () => {
        var ctrGroup = document.createElement("div");
        ctrGroup.className = "wfplot-ctr";

        if(typeof jspdf === "undefined" || typeof saveAs === "undefined") {
          if(!window.navigator.onLine) ctrInernetErr(ctrGroup, el);
          var scripts_loaded = new Set();
          await load_scripts([
            'https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.4.0/jspdf.umd.min.js',
            'https://cdnjs.cloudflare.com/ajax/libs/FileSaver.js/2.0.0/FileSaver.min.js'
          ], scripts_loaded);
          var timeCount = 0;
          await new Promise(resolve => {
            let checkLoaded = setInterval(()=>{
              if(scripts_loaded.size === 2 && timeCount <= 5000) {
                resolve("loaded");
                clearInterval(checkLoaded);
              } else if (timeCount <= 5000){
                timeCount += 500;
              } else {
                clearInterval(checkLoaded);
                ctrInernetErr(ctrGroup, el);
              }
            }, 500);
          })
        }

        ctrGroup.innerHTML =
        `
        <button data-for="png">PNG</button>
        <button data-for="jpg">JPG</button>
        <button data-for="svg_full">SVG (full)</button>
        <button data-for="svg_min">SVG (min)</button>
        <button data-for="pdf">PDF</button>
        <button data-for="dot">Graphviz</button>
        `
        el.appendChild(ctrGroup);
        ctrGroup.querySelector('button[data-for="png"]').addEventListener('click', toPng)
        ctrGroup.querySelector('button[data-for="jpg"]').addEventListener('click', toJpg)
        ctrGroup.querySelector('button[data-for="svg_full"]').addEventListener('click', toSvgFull)
        ctrGroup.querySelector('button[data-for="svg_min"]').addEventListener('click', toSvgMin)
        ctrGroup.querySelector('button[data-for="pdf"]').addEventListener('click', toPdf)
        ctrGroup.querySelector('button[data-for="dot"]').addEventListener('click', ()=>{toDot(dotStr)})
      })();
      return el;
    }

     function addPanZoom(el) {
      (async () => {
        var panError = document.createElement("div");
        panError.className = "panzoom-error";

        if(typeof svgPanZoom === "undefined") {
          if(!window.navigator.onLine) panZoomInernetErr(panError, el);
          var scripts_loaded = new Set();
          await load_scripts([
            'https://cdn.jsdelivr.net/npm/svg-pan-zoom/dist/svg-pan-zoom.min.js'
          ], scripts_loaded);
          var timeCount = 0;
          await new Promise(resolve => {
            let checkLoaded = setInterval(()=>{
              if(scripts_loaded.size === 1 && timeCount <= 5000) {
                resolve("loaded");
                clearInterval(checkLoaded);
              } else if (timeCount <= 5000){
                timeCount += 500;
              } else {
                clearInterval(checkLoaded);
                panZoomInernetErr(panError, el);
              }
            }, 500);
          });

        }
        svgPanZoom(el.querySelector('svg'), {controlIconsEnabled: true});
      })();


      return el;
     }
    /**************/

    return {
      renderValue: function(x) {
        // don't know why Rmd like to append this to table and doesn't give it a </div> closing, remove it
        if (x.rmd) x.dot = x.dot.replaceAll(`<div class='horizontal-scroll'>`, '');
        var plot = document.querySelector(`#${x.plotid}`);

        var viz = new Viz();
        var legendSrc = x.legend_uri === "" ?
          document.querySelector('link[id*="plotwf_legend"]').attributes.href.value :
          x.legend_uri;
        var dotStr = x.dot.replace('plotwf_legend-src\.png', legendSrc);
        [...el.children].forEach(e => e.remove()); //clear all children before attach

        viz[x.plot_method](dotStr, {images: [{path: legendSrc, width: '450px', height: '250px'}]})
        .then(function(plot_el) {
          plot_el.id = x.plotid;
          el.style.width = x.width ? x.width: "100%";
          el.style.height = x.height ? x.height: "100%";
          el.classList.add("wfplot");

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

          if (x.rmd) plot_el.removeAttribute('height width');

          if(x.responsive) makeResponsive(plot_el);
          document.dispatchEvent(new Event('wf_plot_created'));
          return el;
        })
        .then(el=>{
          if(x.plot_ctr) addControl(el, dotStr);
          if(x.pan_zoom) addPanZoom(el);
        })
        .catch(e => {
          var p = document.createElement("pre");
          console.log(e);
          p.innerText = e;
          el.appendChild(p);
        });

      },
      resize: function(width, height) {
        // TODO: code to re-render the widget with a new size
      }
    };
  }
});

/*

*/
