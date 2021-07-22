




HTMLWidgets.widget({
  name: 'plotwf',
  type: 'output',
  factory: function(el, width, height) {
    // TODO: define shared variables for this instance
    return {
      renderValue: function(x) {
        // don't know why Rmd like to append this to table and doesn't give it a </div> closing, remove it
        if (x.rmd) x.dot = x.dot.replaceAll(`<div class='horizontal-scroll'>`, '');
        var plot = document.querySelector(`#${x.plotid}`);
        function makeResponsive(plotid){
         if (plot) {
          if (plot.width) {plot.removeAttribute("width")}
          if (plot.height) {plot.removeAttribute("height")}
          plot.style.width = "100%";
          plot.style.height = "100%";
         }
        }

        var viz = new Viz();
        viz[x.plot_method](x.dot)
        .then(function(plot_el) {
          plot_el.id = x.plotid;
          el.style.width = x.width ? x.width: "100%";
          el.style.height = x.height ? x.height: "100%";
          el.style.overflow = "auto";

          [...el.children].forEach(e => e.remove()); //clear all children before attach
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

          if(x.responsive) makeResponsive(x.plotid);
          document.dispatchEvent(new Event('wf_plot_created'));
        })
        .catch(e => {
          var p = document.createElement("pre");
          console.log(e);
          console.log(x.dot);
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


